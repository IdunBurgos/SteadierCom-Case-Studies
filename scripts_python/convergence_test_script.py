
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from datetime import datetime
from html import escape
from pathlib import Path
import statsmodels.api as sm
import matplotlib.pyplot as plt

rng = np.random.default_rng(42)

import os

import logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(levelname)-8s | %(message)s",
)


def find_samples_statistics_bootstrap(dfs_tot_abs,sample_sizes,cutoff=1e-6):
    """
    Function for making sets of variables and finding SE, and RSE of the samples.
    """

    std = lambda x:np.std(x,ddof=1)
    cv = lambda x: std(x)/ abs(np.mean(x)) * 100 
        
    S = 200  
    batches_mean = {}
    batches_std = {}
    batches_cv = {}
    batches_freq = {}


    population_mean = dfs_tot_abs.mean(axis=1)
    
    population_bool = dfs_tot_abs.values>cutoff
    population_freq = pd.Series(population_bool.mean(axis=1),index=dfs_tot_abs.index)


    n_vars = dfs_tot_abs.shape[0]
    n_cols = dfs_tot_abs.shape[1]
    idx_series = dfs_tot_abs.index

    for n in sample_sizes:
        # Running aggregates across S replicates
        # For replicate means (per variable)

        n = int(n)
        print(n)
        mean_acc = np.zeros(n_vars, dtype=np.float64)
        M2_acc   = np.zeros(n_vars, dtype=np.float64)  # for variance via Welford
        count    = 0

        # For replicate frequencies (per variable): mean of (samples > cutoff)
        freq_acc = np.zeros(n_vars, dtype=np.float64)

        for _ in range(S):
            # sample n columns with replacement (faster + lower overhead than DataFrame.sample)
            col_idx = rng.integers(0, n_cols, size=n)
            # slice as NumPy for speed; rows stay in the same order
            sample_vals = dfs_tot_abs.values[:, col_idx]

            # replicate mean across those n columns (per variable)
            m = sample_vals.mean(axis=1)

            # Welford update for mean and variance of replicate means across S
            count += 1
            delta  = m - mean_acc
            mean_acc += delta / count
            M2_acc   += delta * (m - mean_acc)

            # replicate frequency: fraction > cutoff among the n columns
            f = (sample_vals > cutoff).mean(axis=1)
            freq_acc += f

            # free references ASAP
            del sample_vals

        # finalize across S replicates
        mean_of_means = pd.Series(mean_acc, index=idx_series)
        # unbiased std across replicates of the mean
        if S > 1:
            var_of_means = M2_acc / (S - 1)
        else:
            var_of_means = np.zeros_like(M2_acc)
        std_of_means  = pd.Series(np.sqrt(var_of_means), index=idx_series)

        # CV across replicates of the mean: std / |mean| * 100
        with np.errstate(divide='ignore', invalid='ignore'):
            cv_vals = (std_of_means.values / np.abs(mean_of_means.values)) * 100.0
            cv_vals[~np.isfinite(cv_vals)] = np.nan
        cv_of_means = pd.Series(cv_vals, index=idx_series)

        # average frequency across replicates
        freq_mean = pd.Series(freq_acc / S, index=idx_series)

        batches_mean[n] = mean_of_means
        batches_std[n]  = std_of_means
        batches_cv[n]   = cv_of_means
        batches_freq[n] = freq_mean

    return batches_mean,batches_std,batches_cv,batches_freq,population_mean,population_freq

def plot_point_plot(samples_std_or_cv,ax=None,title="",ylabel=""):

    ax_given = ax is not None

    if ax is None:
        fig,ax = plt.subplots(figsize=(15,5))

    ax = sns.pointplot(pd.DataFrame(samples_std_or_cv),ax=ax, color="blue")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Number of samples")

    if ax_given:
        return ax
    else:
        return fig,ax
    
def plot_mean_vs_cv_plot(datas_mean,datas_cv,population_mean,ax=None,title="",ylim=False,xlim=False,xlabel=""):

    tol = 1e-6
    ax_given = ax is not None
    if ax is None:
        fig,ax = plt.subplots(figsize=(15,5))

    colors = plt.get_cmap("Set3").colors
    i=0

    population_mean_new = population_mean[population_mean>tol]
    for N,mean_mean in datas_mean.items():
        mean_cv = datas_cv[N]

        mean_cv_new = mean_cv[population_mean>tol].map(np.log10)
        mean_mean_new = mean_mean[population_mean>tol].map(np.log10)
        smoothed = sm.nonparametric.lowess(exog=mean_mean_new, endog=mean_cv_new, frac=0.2)  
        
        ax.plot(10**smoothed[:, 0], 10**smoothed[:, 1], c=colors[i],label = str(N))
        i+=1

    hist, bins = np.histogram(population_mean_new.values, bins=10)
    
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    ax_twin = ax.twinx()
    ax_twin.hist(population_mean_new.values,bins=logbins,alpha=0.1,label="Number of variables")
    ax_twin.set_ylabel("Number of variables")
    ax_twin.set_yscale("log")
    ax_twin.legend(title="Histogram of total 'population'", loc="upper left")

    ax.set_xlabel(xlabel=xlabel)
    ax.set_ylabel("Uncertainty in mean (RSE %)")
    ax.set_title(title)
    ax.legend()
    ax.set_xscale("log")
    
    if ylim:
        ax.set_ylim(ylim)
    if xlim:
        ax.set_xlim(xlim)
    
    if ax_given:
        return ax
    else:
        return fig,ax

def highlight_row(x):
    color = 'background-color: yellow'
    df1 = pd.DataFrame('', index=x.index, columns=x.columns)
    row1_index = x[x["Change in RSE [%]"]<-1].index[-1]
    row2_index = x[x["Change in RSE [%]"]>-1].index[0]

    df1.loc[row1_index, :] = color
    df1.loc[row2_index, :] = color
    return df1

def evaluate_improvement_sample_size(datas_std_or_cv,population_mean,population_freq):
    stats = pd.DataFrame.from_dict(datas_std_or_cv)
    stats = stats[(population_mean>population_mean.mean()*1e-2) & (population_freq>1e-2)].mean()

    diff = pd.DataFrame(stats.T.iloc[1:].reset_index(drop=True)-stats.T.iloc[:-1].reset_index(drop=True))
    sample_sizes = list(datas_std_or_cv.keys())
    index_a = sample_sizes[:-1]
    index_b = sample_sizes[1:]
    index = pd.Index(list(zip(index_a,index_b)))
    diff.set_index(keys=index,inplace=True)
    diff_mean = diff.T.mean()
    diff_mean.name="Change in RSE [%]"

    current_rse = stats.T.iloc[1:]
    current_rse.index=index
    current_rse.name = "Current RSE [%]"
    combined = pd.concat([diff_mean,current_rse],axis=1)
    combined.style.apply(highlight_row,axis=None)
    return combined

def load_data_and_preprocess(output_dir):
    logging.info(f"Loading data from {output_dir}...")
    dfs = []
    for d in os.listdir(output_dir):
        if d.endswith(".tsv"):
            df = pd.read_csv(output_dir + str(d),sep="\t",index_col=(0,1))
            dfs.append(df)

    dfs_tot = pd.concat(dfs,axis=1)
    dfs_tot = dfs_tot.T.drop_duplicates().T
    dfs_tot.columns = np.arange(0,dfs_tot.shape[1])

    # Separate positive and degative values
    pos = dfs_tot[dfs_tot>0].dropna(how="all").copy()
    neg = dfs_tot[dfs_tot<0].dropna(how="all").map(abs).copy()

    neg["type"]="neg"
    pos["type"]="pos"
    neg.set_index("type",append=True,inplace=True)
    pos.set_index("type",append=True,inplace=True)

    dfs_tot_abs_nan = pd.concat([pos,neg],axis=0)
    dfs_tot_abs = dfs_tot_abs_nan.fillna(0)

    return dfs_tot, dfs_tot_abs

def process_data(dfs_tot_abs):

    logging.info(f"Processing data (finding standard deviation and coefficient of variance)...")

    sample_sizes = np.array([250*(2**i) for i in np.arange(0,6)])
    sample_sizes_higher_res = np.arange(500,5500,500)
    
    datas_mean,datas_std,datas_cv,datas_freq,population_mean,population_freq = find_samples_statistics_bootstrap(dfs_tot_abs,sample_sizes)
    datas_high_res_mean,datas_high_res_std,datas_high_res_cv,datas_freq_high_res,population_mean,population_freq = find_samples_statistics_bootstrap(dfs_tot_abs,sample_sizes_higher_res)
    
    results = {"dfs_tot_abs":dfs_tot_abs,
               "mean_std_cv":(datas_mean,datas_std,datas_cv,datas_freq),
               "mean_std_cv_high_res":(datas_high_res_mean,datas_high_res_std,datas_high_res_cv,datas_freq_high_res),
               "population_mean": population_mean,
               "population_freq": population_freq
               }
    return results

def make_report(results,output_dir,ecoli_or_gut="Venturelli2018"):

    html_dir = os.path.join("html_reports/",ecoli_or_gut,output_dir.split("/")[-2])
    fig_dir = os.path.join(html_dir, "figures")

    os.makedirs(fig_dir,exist_ok=True)

    dfs_tot_abs = results["dfs_tot_abs"]
    datas_mean,datas_std,datas_cv,datas_freq = results["mean_std_cv"]
    datas_high_res_mean,datas_high_res_std,datas_high_res_cv,datas_freq_high_res = results["mean_std_cv_high_res"]
    population_mean = results["population_mean"]
    population_freq = results["population_freq"]

    mean_of_mean = dfs_tot_abs.T.mean().mean()
    max_of_mean = dfs_tot_abs.T.mean().max()
    median_of_mean = dfs_tot_abs.T.mean().median()

    overall_results = pd.Series({"Mean of 'population' means":mean_of_mean,
                "Max of 'population' means":max_of_mean,
                "Median of 'population' means":median_of_mean               
               })
    
    overall_results_df = pd.DataFrame(overall_results) 
    overall_results_df.columns = ["value"]
    overall_results_df.index.names = ["metric"]
    
    improvement_df = evaluate_improvement_sample_size(datas_high_res_cv,population_mean,population_freq)
    improvement_df.index.names = ["previoius sample size","current sample size"]
   
    fig_std_cv,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    ax1 = plot_point_plot(datas_std,
                         ax=ax1, 
                         title = "SE at various sample sizes",
                         ylabel="SE")
    ax2 = plot_point_plot(datas_cv, 
                         ax=ax2,
                         title = "RSE at various sample sizes",
                         ylabel="RSE [%]")

    fig_std_cv.tight_layout()
    fig_std_cv.savefig(os.path.join(fig_dir,"fig_std_cv.png"))
    plt.close(fig_std_cv)

    fig_flux_sample,(ax1,ax2) = plt.subplots(1,2,figsize=(14,5))
    ax1 = plot_mean_vs_cv_plot(datas_mean,datas_cv,population_mean,ax=ax1,title="Mean rate vs RSE",xlabel="Mean value",xlim=[1e-7,10])
    ax2 = plot_mean_vs_cv_plot(datas_freq,datas_cv,population_freq,ax=ax2,title="Frequency vs RSE",xlabel="Frequency of observation")

    fig_flux_sample.tight_layout()
    fig_flux_sample.savefig(os.path.join(fig_dir,"fig_flux_sample.png"))
    plt.close(fig_flux_sample)
    
    
    r = HTMLReport("Convergence analysis","Idun Burgos")

    r.add_heading(f"General results across {dfs_tot_abs.shape[1]} samples")
    r.add_paragraph("As a basis for this analysis we calculate some general properties of our data. \n"\
                    f"We have <b>{dfs_tot_abs.shape[1]} samples</b> and <b>{dfs_tot_abs.shape[0]} metabolite-member interactions (variables)</b>, which is assumed to be a larger number than what in most cases will be necessary to perform." \
                    " This therefore serves as an approximate to the 'population'.\n" \
                    "<b>SteadierCom</b> takes a number of samples, a batch, and gives the user the mean value of each variable from this batch. Based on the sample mean of "\
                    "the 'population' dataset described above, the following results show the max, median and mean flux of the variables:"
                    )
    r.add_table(overall_results_df)

    r.add_heading("Bootstrapping for calculating SE and RSE")
    r.add_paragraph("In order to evaluate the effect of increasing sample size we have created <i><b>S</b></i> batches of <i><b>N</b></i> samples. We assess the similarity between these batches "\
                    "by computing the <b>standard error (SE)</b> of the batch means for each variable. Since the mean of each variable is quite different we also look at "\
                    "the <b>realtive SE (RSE)</b>, which expresses the variability of the mean of each variable relative to their mean."
                    )
    
    r.add_paragraph("The <i><b>S</b></i>=200 batches were created using bootstrapping with replacement at various <i><b>N</b></i> values."
                    )
    r.add_heading("SE and RSE vs batch size",level=3)    

    r.add_paragraph("This figure shows the change in standard error of the mean (SE) and relative standard error (RSE) when batch size increases (has more samples)."\
                    " A low SE means that the results from simulations with a similar number of samples are quite similar.")

    r.add_image(os.path.join("figures",
                             "fig_std_cv.png"), 
                            caption="Change in SE at various batch sizes", 
                            width=None)

    r.add_heading("RSE at different means and frequencies of variable observations",level=3)
    r.add_paragraph("As described previously, the metabolite-member interactions (variables) will have different mean values. It is also expected that "\
                    "some samples will give a non-zero value of a metabolite-member interactions, while others will not. From this we can both calculate the mean of the "\
                    "variable in the batch and also the frequency that the interaction was observed. In this part of the report, study RSE as a function of the mean flux across"\
                    " all batches (mean of means, left plot) and as a function of frequency across all batches (mean of frequencies, right plot)."
                    )
    
    r.add_paragraph("As described previously, the metabolite-member interactions (variables) will have different mean values. It is also expected that "\
                    "some samples will give a non-zero value of a metabolite-member interactions, while others will not. From this we can both calculate the mean of the "\
                    "variable in the batch and also the frequency that the interaction was observed. In this part of the report, study RSE as a function of the mean flux across"\
                    " all batches (mean of means, left plot) and as a function of frequency across all batches (mean of frequencies, right plot)."
                    )

    r.add_image(os.path.join("figures",
                             "fig_flux_sample.png"), 
                            caption="Change in RSE at various sample sizes. The lines are made with statsmodel lowess and represent the scatterplot created by plotting the "\
                                    f" mean and RSE of each variable at different batch sizes. The histogram shows the underlying distribution in the {dfs_tot_abs.shape[1]} "\
                                    "samples, i.e. ~population mean or frequency of each variable.", 
                            width=None)
    
    r.add_heading("Effect of increasing sample size by 500 on RSE [%]")
    r.add_paragraph("From the previous plot we see that how the frequency of an observed interaction (variable) and its mean value affect the RSE value. "\
                    "The results show that there is more relative uncertainity in the interactions that rarely occur or occur with a low mean rate. Moreover, "\
                    " frequency and mean are related, since lower frequency of an interaction also drives down its mean value.")

    r.add_paragraph("To account for this in determining sample size, we will filter values based on ~population median and frequency. Interactions (variables) that "\
                    "occur with <1% of samples (<0.01 frequency) in the population are unlikely to be of importance. The same is true for interactions that occur at a"\
                    " very low mean rate for example <1% of the mean of the population means of all variables.")
    r.add_table(improvement_df)

    r.save(os.path.join(html_dir,"report.html"))



class HTMLReport:
    def __init__(self, title="Report", author=None):
        self.title = title
        self.author = author or ""
        self.blocks = []


    def _wrap(self):
        return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{escape(self.title)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  :root {{ --fg:#111; --muted:#666; --bg:#fff0f5; --card:#fafafa; --accent:#0b5fff; }}
  .report {{
    font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
    color: var(--fg);
    background: var(--bg);
    margin: 24px;
  }}
  header {{ margin-bottom: 24px; }}
  h1 {{ margin: 0 0 4px 0; }}
  h2 {{ margin-top: 28px; }}
  p, li {{ line-height: 1.55; }}
  .meta {{ color: var(--muted); font-size: 0.9rem; }}
  .card {{ background:var(--card); border:1px solid #eee; border-radius:12px; padding:16px; margin: 14px 0; }}
  .img {{ text-align:center; }}
  .img img {{ max-width:100%; height:auto; border-radius:12px; }}
  figure {{ margin:0; }}
  figcaption {{ color:var(--muted); font-size:0.9rem; margin-top:6px; }}
  table {{ border-collapse: collapse; width:50%; }}
  th, td {{ border:1px solid #e6e6e6; padding:8px 10px; text-align:left; }}
  th {{ background:#f4f6ff; }}
  code {{ background:#f6f8fa; padding:2px 6px; border-radius:6px; }}
  .kpi {{ display:inline-block; padding:12px 14px; margin:6px 8px 6px 0; background:#f7f7f7; border:1px solid #eee; border-radius:12px; }}
  .kpi b {{ font-size:1.15rem; }}
</style>
</head>
<body>
  <div class="report">
    <header>
      <h1>{escape(self.title)}</h1>
      <div class="meta">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} {("· " + escape(self.author)) if self.author else ""}</div>
    </header>
    {"".join(self.blocks)}
  </div>
</body>
</html>"""

    def add_heading(self, text, level=2):
        self.blocks.append(f"<h{level}>{escape(text)}</h{level}>")

    def add_paragraph_old(self, text):
        self.blocks.append(f"<p>{escape(text)}</p>")

    def add_paragraph(self, text, escape_html=True):
        if escape_html:
            from html import escape
            text = escape(text)
        self.blocks.append(f"<p>{text}</p>")

    def add_kpi(self, label, value):
        self.blocks.append(f'<div class="kpi"><div>{escape(label)}</div><b>{escape(str(value))}</b></div>')

    def add_image(self, path, caption=None, width=None):
        path = str(path)
        style = f' style="max-width:{int(width)}px"' if width else ""
        cap = f"<figcaption>{escape(caption)}</figcaption>" if caption else ""
        self.blocks.append(f'''
        <div class="card img">
          <figure>
            <img src="{escape(path)}"{style} loading="lazy" />
            {cap}
          </figure>
        </div>
        ''')

    def add_table(self, data, caption=None):
        # data: pandas.DataFrame OR list[dict] OR list[list] with header
        html = None
        try:
            if isinstance(data, pd.io.formats.style.Styler):
                # styled DataFrame
                html = data.to_html()
            elif hasattr(data, "to_html"):
                html = data.to_html(border=0, classes="table", index=True)
        except Exception:
            pass

        if html is None:
            # Fallbacks for lists/dicts
            rows = []
            if isinstance(data, list) and data and isinstance(data[0], dict):
                headers = list({k for row in data for k in row.keys()})
                head = "".join(f"<th>{escape(str(h))}</th>" for h in headers)
                for row in data:
                    cells = "".join(f"<td>{escape(str(row.get(h,'')))}</td>" for h in headers)
                    rows.append(f"<tr>{cells}</tr>")
            elif isinstance(data, list) and data and isinstance(data[0], (list, tuple)):
                headers = data[0]
                head = "".join(f"<th>{escape(str(h))}</th>" for h in headers)
                for row in data[1:]:
                    cells = "".join(f"<td>{escape(str(c))}</td>" for c in row)
                    rows.append(f"<tr>{cells}</tr>")
            else:
                # last resort: single-column
                head = "<th>Value</th>"
                rows = [f"<tr><td>{escape(str(x))}</td></tr>" for x in (data if isinstance(data, list) else [data])]
            html = f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"

        card = f'<div class="card">{html}'
        if caption:
            card += f"<div class='meta' style='margin-top:8px'>{escape(caption)}</div>"
        card += "</div>"
        self.blocks.append(card)

    def add_html(self, raw_html):
        self.blocks.append(raw_html)

    def save(self, path="report.html"):  
        Path(path).write_text(self._wrap(), encoding="utf-8")
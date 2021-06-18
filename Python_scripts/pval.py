import argparse
import logging.config
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats._continuous_distns import _distn_names  # import distributions names
import warnings
import matplotlib
import matplotlib.pyplot as plt
import re
from statsmodels.sandbox.stats.multicomp import multipletests
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

parser = argparse.ArgumentParser()
statis = importr('stats')

parser.add_argument('-df', '--dataframe', help= 'Fold Change dataframe')
parser.add_argument('-o', '--output', help= 'Output directory' )
parser.add_argument('-n', '--name', help= 'file name')
parser.add_argument('-a', '--all_in_one', help = 'measure pvalues building one big distribution with all Fold changes',
                    action= "store_true")
parser.add_argument('-p', '--per_row', help = 'measure pvalues building one distribution per row', action='store_true')
parser.add_argument('-c', '--comparing', help = 'measure pvalues per comparisons', action= 'store_true')
parser.add_argument('-s', '--significant', help = 'path + filename to obtain significant pvalues element')
parser.add_argument('-k', '--correction', help= 'get a pvalue method correction (for -p and -c only) through the '
                                                 'following list : '
                                                 '"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", '
                                                 '"fdr", "none"')

args = parser.parse_args()




def compute_z_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add one column with z-score of ratio already computed
    """
    # TODO: if no columns named "ratio" / if several columns with "ratio_"
    df['zscore'] = stats.zscore(df)
    return df


def find_best_distribution(df: pd.DataFrame, out_histogram_distribution: str):
    """
    Find best distribution among all the scipy.stats distribution and returns it with its parameters
    """
    dist = np.around(np.array((df['zscore']).astype(float)), 2)  # TODO : pourquoi arrondir ??

    best_dist, best_dist_name, best_fit_params = get_best_fit(dist, out_histogram_distribution)
    print(best_dist_name)
    print(best_fit_params)

    # logger.info("Best fit is", str(best_dist_name), "with", str(*best_fit_params))
    args_param = dict(e.split('=') for e in best_fit_params.split(', '))
    for k, v in args_param.items():
        args_param[k] = float(v)

    best_distribution = getattr(stats, best_dist_name)
    q_val = best_dist.ppf(0.95, **args_param)
    print("And the q value is", q_val)
    return best_distribution, args_param


def get_best_fit(input_array, out_file):
    matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
    matplotlib.style.use('ggplot')
    """Return the best fit distribution to data and its parameters"""

    # Load data
    data = pd.Series(input_array)

    # Find best fit distribution
    best_fit_name, best_fit_params = best_fit_distribution(data, 200)

    best_dist = getattr(stats, best_fit_name)

    # Make probability density function (PDF) with best params
    pdf = make_pdf(best_dist, best_fit_params)

    # parameters
    param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
    print(best_fit_params)

    # get significant numbers to keep for distribution parameters
    sn = 3

    # !!!!! Attention si la distribution a un paramètres qui ne peut pas être égal à 0 par définition mais se retrouve avec cette valeur après le .format()
    param_str = ', '.join(['{}={:0.{}f}'.format(k, v, sn) for k, v in zip(param_names, best_fit_params)])
    dist_str = '{} ({})'.format(best_fit_name, param_str)

    print("dist_str : ", dist_str)

    # # Display
    param_str_plot = ', '.join(['{}={:0.2f}'.format(k, v) for k, v in zip(param_names, best_fit_params)])
    dist_str_plot = '{} ({})'.format(best_fit_name, param_str_plot)
    plot_best_fit(data, dist_str_plot, pdf, out_file)

    return best_dist, best_fit_name, param_str


def plot_best_fit(data, dist_str, pdf, out_file):
    # #plot methode alternative https://stackoverflow.com/questions/6620471/fitting-empirical-distribution-to-theoretical-ones-with-scipy-python

    plt.figure(figsize=(12, 8))
    plt.hist(data, bins=150, density=True, alpha=0.5, label='Data')
    plt.plot(pdf, lw=2, label='PDF')
    plt.legend(loc='upper right', shadow=True, fontsize='x-large')
    plt.title(u'Best fit distribution \n' + dist_str)
    plt.xlabel(u'z-score')
    plt.ylabel('frequency')
    plt.savefig(out_file)
    return


def make_pdf(dist, params, size=10000):
    """Generate distributions's Probability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf


def get_scipy_distributions():
    distribution_list = []
    print(_distn_names)
    for dist in _distn_names:
        distribution_list.append(dist)

        # distributions toutes assez rapides à calculer, pour vue rapide des résultats dégager "levy_stable"

    stats_distributions = [getattr(stats, d) for d in distribution_list]

    return stats_distributions


def best_fit_distribution(data, bins=200):
    notlevy = re.compile(r'.*levy_stable.*')
    matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
    matplotlib.style.use('ggplot')

    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)

    x = (x + np.roll(x, -1))[:-1] / 2.0

    DISTRIBUTIONS = get_scipy_distributions()
    print(len(DISTRIBUTIONS))
    # Best holders
    best_distribution = stats.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in DISTRIBUTIONS:
        print(distribution)
        name = str(distribution)
        if re.match(notlevy, name) is None:
            # Try to fit the distribution
            try:
                # Ignore warnings from data that can't be fit
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')

                    # fit dist to data
                    params = distribution.fit(data)

                    # Separate parts of parameters
                    arg = params[:-2]
                    loc = params[-2]
                    scale = params[-1]

                    # Calculate fitted PDF and error with fit in distribution
                    pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                    sse = np.sum(np.power(y - pdf, 2.0))
                    # identify if this distribution is better
                    if best_sse > sse > 0:
                        best_distribution = distribution
                        best_params = params
                        best_sse = sse
                        print("Distribution : {}" / format(distribution))
                        print(arg)

            except Exception:
                pass

    return best_distribution.name, best_params


def compute_p_value(df: pd.DataFrame, test: str, best_dist, args_param) -> pd.DataFrame:
    """
    compute p-value depending on the chosen test
    """
    if test == 'right-tailed':
        df['pvalue'] = 1 - best_dist.cdf(df['zscore'], **args_param)
    elif test == 'left-tailed':
        df['pvalue'] = best_dist.cdf(df['zscore'], **args_param)
    elif test == 'two-sided':
        df['pvalue_right'] = 1 - best_dist.cdf(df['zscore'], **args_param)
        df['pvalue_left'] = best_dist.cdf(df['zscore'], **args_param)

        for i in df.index.values:
            pr = df.loc[i, 'pvalue_right']
            pl = df.loc[i, 'pvalue_left']
            df.loc[i, 'pvalue'] = 2 * (min(pr, pl))
            df.loc[i, 'pvalue_best_side'] = min(pr, pl)

    else:
        print("WARNING: two-tailed / left-tailed / right-tailed ")
    return df


def pvalues_all_in_one(df, out):
    output = pd.DataFrame(index=df.index)
    longlist = list()
    n = len(df.index)
    for i, c in enumerate(df):
        dtmp = df[c]
        for z in dtmp.values:
            longlist.append(z)
    inter = longlist.copy()
    allin = pd.DataFrame()
    outlist = list()
    for z, e in enumerate(longlist):
        if e > 10 or e < 1E-4:
            outlist.append(z)
    outlist.sort(reverse= True)
    for o in outlist:
        del inter[o]
    for o in outlist:
        longlist[o] = np.mean(inter)

    allin['FC'] = longlist

    tmp = compute_z_score(pd.DataFrame(allin))
    fit, param = find_best_distribution(tmp, args.output)
    pval = compute_p_value(tmp, 'two-sided', fit, param)
    mat = pval['pvalue_best_side'].values.tolist()
    for z in outlist:
        mat[z] = 0
    for i, c in enumerate(df):
        output[c] = mat[0:n]
        del mat[0:n]


    output.to_csv(out)
    return output


def pvalues_per_comparing(df, out):
    output = pd.DataFrame(index=df.index)
    for i, c in enumerate(df):
        print(i)
        dtmp = pd.DataFrame(df[c])
        tmp = compute_z_score(dtmp)
        fit, param = find_best_distribution(tmp, args.output)
        pval = compute_p_value(tmp, 'two-sided', fit, param)

        if args.correction is not None:
            pval['pvalue_best_side'] = statis.p_adjust(FloatVector(pval['pvalue_best_side']), method = args.correction)
            output[c] = pval['pvalue_best_side']
        else:
            output[c] = pval['pvalue_best_side']

    output.to_csv(out)
    return output


def pvalues_per_rows(df, out):
    output = pd.DataFrame(index=df.index)
    for i, r in enumerate(df.index):
        print(i)
        outlist = list()
        dtmp = pd.DataFrame(df.loc[r])
        for z, e in enumerate(dtmp[r].tolist()) :
            if e > 300 or e < 1E-4:
                outlist.append(dtmp.index[z])
                dtmp = dtmp.drop(dtmp.index[z])

        tmp = compute_z_score(dtmp)
        fit, param = find_best_distribution(tmp, args.output)
        pval = compute_p_value(tmp, 'two-sided', fit, param)
        if i == 0:
            for u in pval.index :
                output[u] = np.NaN
        for d in outlist:
            newrow = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
            pval.loc[d] = newrow
        if args.correction is not None:
            pval['pvalue_best_side'] = statis.p_adjust(FloatVector(pval['pvalue_best_side']), method = args.correction)
        output.loc[r] = pval['pvalue_best_side']

    output.to_csv(out)
    return output


def sig_exctract(pvalues:pd.DataFrame):
    sigs = pd.DataFrame()
    for c in pvalues:
        val_low = pvalues[pvalues[c] < 0.05][c]

        sigs = pd.concat([sigs, pd.Series(val_low)], ignore_index=True, axis = 1)
    sigs.columns = pvalues.columns

    sigs.to_csv(args.significant)


if __name__ == "__main__":
    df = pd.read_csv(args.dataframe, index_col= 0)
    out = args.output + args.name

    if args.all_in_one:
        tmp = pvalues_all_in_one(df, out)
    if args.per_row:
        tmp =pvalues_per_rows(df, out)
    if args.comparing:
        tmp = pvalues_per_comparing(df, out)

    if args.significant is not None:
        sig_exctract(tmp)





    # dat = pd.read_csv("/run/media/aurelien/ACOFFE/Stage/essais/all_FC.csv", index_col=0 )
    #
    # output = pd.DataFrame(index= dat.index)
    # for i, c in enumerate(dat):
    #     print(i)
    #     dtmp = pd.DataFrame(dat[c])
    #     tmp = compute_z_score(dtmp)
    #     fit, param = find_best_distribution(tmp, "/run/media/aurelien/ACOFFE/Stage/essais/")
    #     pval = compute_p_value(tmp, 'two-sided', fit, param)
    #     output[c] = pval['pvalue_best_side']
    # output.to_csv("/run/media/aurelien/ACOFFE/Stage/essais/pvalues.csv")

    # uwu = re.compile(r'T24_AB')
    # for i, c in enumerate(dat.columns):
    #     if re.match(uwu, c):
    #         dtmp = pd.DataFrame(dat[c])
    #         tmp = compute_z_score(dtmp)
    #         fit, param = find_best_distribution(tmp, "/run/media/aurelien/ACOFFE/Stage/essais/")
    #         pval = compute_p_value(tmp, 'two-sided', fit, param)

    # dat = pd.read_csv("/run/media/aurelien/ACOFFE/Stage/essais/all_FC_tot.csv", index_col=0 )
    # output = pd.DataFrame(index= dat.index)
    # longlist = list()
    # n = len(dat.index)
    # for i, c in enumerate(dat):
    #     dtmp = dat[c]
    #     for i in dtmp.values :
    #         longlist.append(i)
    #
    # allin = pd.DataFrame()
    # allin['FC'] = longlist
    #
    # tmp = compute_z_score(pd.DataFrame(allin))
    # fit, param = find_best_distribution(tmp, "/run/media/aurelien/ACOFFE/Stage/essais/")
    # pval = compute_p_value(tmp, 'two-sided', fit, param)
    # mat = pval['pvalue_best_side'].values.tolist()
    # for i, c in enumerate(dat):
    #     output[c] = mat[0:n]
    #     del mat[0:n]
    # print(output)
    # output.to_csv("/run/media/aurelien/ACOFFE/Stage/essais/pvalues_totlog2.csv")


    #     if i == 0:
    #         for u in pval.index :
    #             output[u] = np.NaN
    #     output.loc[c] = pval['pvalue_best_side']
    #     print(output)
    #
    # output.to_csv("/run/media/aurelien/ACOFFE/Stage/essais/pvalues_V2.csv")




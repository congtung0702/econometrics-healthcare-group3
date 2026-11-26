"""
Assignment 6 - Risk prediction (complete analysis, questions A-H)

Put this script in the same folder as Fake_CRC_riskdata.csv and run it.
Needed packages (install once in the PyCharm terminal):
    pip install pandas numpy matplotlib statsmodels scikit-learn

Output:
    - output_risk_prediction/output.txt   (all numbers for Overleaf)
    - output_risk_prediction/*.png        (all figures)
"""
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm
from sklearn.metrics import roc_auc_score, roc_curve

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
FOLDER = Path(__file__).parent
OUT = FOLDER / "output_risk_prediction"
OUT.mkdir(exist_ok=True)
REPORT = open(OUT / "output.txt", "w", encoding="utf-8")

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 20)

SEED = 2026


def out(*args):
    """Print to the console and to output.txt."""
    text = " ".join(str(a) for a in args)
    print(text)
    REPORT.write(text + "\n")


def header(title):
    out("\n" + "=" * 78)
    out(title)
    out("=" * 78)


data = pd.read_csv(FOLDER / "Fake_CRC_riskdata.csv")

# Transformed FIT variables (used in H)
data["log_FIT_1"] = np.log1p(data["FIT_1"])   # log(1 + FIT), defined for FIT = 0
data["log_FIT_2"] = np.log1p(data["FIT_2"])

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def fit_logit(formula, df):
    """Logistic regression (same as glm(..., family = binomial) in R)."""
    return smf.glm(formula, data=df, family=sm.families.Binomial()).fit()


def coefficient_table(model):
    table = pd.DataFrame({
        "Estimate": model.params,
        "Std. Error": model.bse,
        "z-value": model.tvalues,
        "p-value": model.pvalues,
        "Odds Ratio": np.exp(model.params),
    })
    return table


def calibration(actual, predicted, name, filename, bins=10):
    """Calibration plot (10 equally sized bins), intercept and slope with Wald tests."""
    df = pd.DataFrame({"actual": actual, "predicted": predicted})
    df["bin"] = pd.qcut(df["predicted"].rank(method="first"), bins, labels=False) + 1
    table = df.groupby("bin").agg(
        predicted=("predicted", "mean"),
        observed=("actual", "mean"),
        n=("actual", "size"),
    )

    # Plot
    max_value = max(table["predicted"].max(), table["observed"].max()) * 1.05
    plt.figure(figsize=(6, 5))
    plt.scatter(table["predicted"], table["observed"], color="C0", zorder=3)
    plt.plot([0, max_value], [0, max_value], "r--", label="Perfect calibration")
    plt.xlim(0, max_value)
    plt.ylim(0, max_value)
    plt.xlabel("Predicted CRC probability")
    plt.ylabel("Observed CRC proportion")
    plt.title(f"Calibration plot - {name}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / filename, dpi=300)
    plt.close()

    # Calibration intercept and slope: regress outcome on the linear predictor
    p = np.clip(predicted, 1e-8, 1 - 1e-8)
    lp = np.log(p / (1 - p))
    cal_model = sm.GLM(actual, sm.add_constant(lp), family=sm.families.Binomial()).fit()
    est = np.asarray(cal_model.params)
    se = np.asarray(cal_model.bse)
    h0 = np.array([0.0, 1.0])
    z = (est - h0) / se
    pval = 2 * norm.sf(np.abs(z))
    results = pd.DataFrame(
        {"Estimate": est, "Std. Error": se, "H0": h0, "z-value": z, "p-value": pval},
        index=["Calibration intercept", "Calibration slope"],
    )
    return table, results


def roc_plot(curves, filename, title):
    """curves: list of (name, actual, predicted). Returns dict with AUROCs."""
    aucs = {}
    plt.figure(figsize=(6, 5))
    for name, actual, predicted in curves:
        fpr, tpr, _ = roc_curve(actual, predicted)
        auc = roc_auc_score(actual, predicted)
        aucs[name] = auc
        plt.plot(fpr, tpr, lw=2, label=f"{name} (AUROC = {auc:.4f})")
    plt.plot([0, 1], [0, 1], "r--", lw=1)
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.title(title)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(OUT / filename, dpi=300)
    plt.close()
    return aucs


# ---------------------------------------------------------------------------
# A. Descriptive statistics
# ---------------------------------------------------------------------------
header("A. DESCRIPTIVE STATISTICS")
vars_ = ["CRC", "age", "male", "FIT_1", "FIT_2"]
out("Number of observations:", len(data))
out("Missing values per variable:\n" + data[vars_].isna().sum().to_string())

desc = pd.DataFrame({
    "Mean": data[vars_].mean(),
    "SD": data[vars_].std(),
    "Median": data[vars_].median(),
    "Min": data[vars_].min(),
    "Max": data[vars_].max(),
})
out("\nDescriptive statistics:\n" + desc.round(3).to_string())

out(f"\nCRC cases: {int(data['CRC'].sum())} ({data['CRC'].mean() * 100:.2f}%)")
out(f"Share FIT_1 = 0: {(data['FIT_1'] == 0).mean() * 100:.1f}%")
out(f"Share FIT_2 = 0: {(data['FIT_2'] == 0).mean() * 100:.1f}%")
both_zero = (data["FIT_1"] == 0) & (data["FIT_2"] == 0)
out(f"Share both FIT = 0: {both_zero.mean() * 100:.1f}%  "
    f"(CRC rate {data.loc[both_zero, 'CRC'].mean() * 100:.2f}% vs "
    f"{data.loc[~both_zero, 'CRC'].mean() * 100:.2f}% for the others)")

out("\nAge distribution and CRC rate per age:")
age_table = data.groupby("age")["CRC"].agg(N="size", CRC_rate="mean")
age_table["CRC_rate"] = (age_table["CRC_rate"] * 100).round(2)
out(age_table.to_string())

out("\nCRC rate per sex (%):")
out((data.groupby("male")["CRC"].mean() * 100).round(2).to_string())

out("\nCRC rate (%) by previous FIT category:")
bins = [-1, 0, 10, 20, 30, 47]
labels = ["0", "0-10", "10-20", "20-30", "30-47"]
fit_table = pd.DataFrame()
for fit in ["FIT_1", "FIT_2"]:
    g = data.groupby(pd.cut(data[fit], bins, labels=labels), observed=True)["CRC"]
    fit_table[f"{fit} N"] = g.size()
    fit_table[f"{fit} CRC (%)"] = (g.mean() * 100).round(2)
out(fit_table.to_string())

out("\nCorrelation matrix:")
out(data[vars_].corr().round(3).to_string())

# Figure: CRC rate by FIT category
fig, ax = plt.subplots(figsize=(7, 4.5))
x = np.arange(len(labels))
ax.bar(x - 0.2, fit_table["FIT_1 CRC (%)"], 0.4, label="FIT 1")
ax.bar(x + 0.2, fit_table["FIT_2 CRC (%)"], 0.4, label="FIT 2")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_xlabel("Previous FIT result (µg Hb/g stool)")
ax.set_ylabel("CRC detected in round 3 (%)")
ax.set_title("Observed CRC rate by previous FIT result")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "A_crc_rate_by_fit.png", dpi=300)
plt.close(fig)

# ---------------------------------------------------------------------------
# B & C. Train-test split and baseline logistic regression
# ---------------------------------------------------------------------------
header("B & C. BASELINE LOGISTIC REGRESSION (70% training data)")
rng = np.random.default_rng(SEED)
train_idx = rng.choice(len(data), size=int(0.7 * len(data)), replace=False)
train = data.iloc[train_idx].copy()
test = data.drop(index=data.index[train_idx]).copy()
out(f"Training set: {len(train)} observations ({int(train['CRC'].sum())} CRC cases)")
out(f"Test set:     {len(test)} observations ({int(test['CRC'].sum())} CRC cases)")

f_base = "CRC ~ male + age + FIT_1 + FIT_2"
m_base = fit_logit(f_base, train)
out("\n" + coefficient_table(m_base).round(4).to_string())
out(f"\nOdds ratio per 10 µg/g increase: FIT_1 = {np.exp(10 * m_base.params['FIT_1']):.3f}, "
    f"FIT_2 = {np.exp(10 * m_base.params['FIT_2']):.3f}")
out(f"95% CI for age coefficient: [{m_base.conf_int().loc['age', 0]:.4f}, "
    f"{m_base.conf_int().loc['age', 1]:.4f}]")
out(f"Predicted probability for the intercept (female, age 0, FIT = 0): "
    f"{1 / (1 + np.exp(-m_base.params['Intercept'])) * 100:.4f}%")

test["p_base"] = m_base.predict(test)

# ---------------------------------------------------------------------------
# D. Calibration
# ---------------------------------------------------------------------------
header("D. CALIBRATION (30% test data) - BASELINE")
cal_tab, cal_res = calibration(test["CRC"].values, test["p_base"].values,
                               "Baseline model", "D_calibration_baseline.png")
cal_tab_print = cal_tab.copy()
cal_tab_print[["predicted", "observed"]] = (cal_tab_print[["predicted", "observed"]] * 100).round(3)
out("Calibration per bin (predicted and observed in %):\n" + cal_tab_print.to_string())
out("\n" + cal_res.round(4).to_string())

# ---------------------------------------------------------------------------
# E. ROC and AUROC
# ---------------------------------------------------------------------------
header("E. ROC CURVE AND AUROC - BASELINE")
auc_base = roc_plot([("Baseline", test["CRC"], test["p_base"])],
                    "E_roc_baseline.png", "ROC curve - Baseline model")["Baseline"]
out(f"AUROC baseline model: {auc_base:.4f}")

# ---------------------------------------------------------------------------
# H. Improved models
# ---------------------------------------------------------------------------
header("H. IMPROVED MODELS")
formulas = {
    "Baseline": f_base,
    "Interactions": "CRC ~ (male + age + FIT_1 + FIT_2)**2",
    "log(1+FIT)": "CRC ~ male + age + log_FIT_1 + log_FIT_2",
}

models = {}
for name, f in formulas.items():
    models[name] = fit_logit(f, train)
    test[f"p_{name}"] = models[name].predict(test)

out("Coefficients of the log(1+FIT) model:")
out(coefficient_table(models["log(1+FIT)"]).round(4).to_string())

out("\nCoefficients of the interaction model:")
out(coefficient_table(models["Interactions"]).round(4).to_string())

# Calibration plots and comparison on the same split
rows = []
for name in formulas:
    fname = "H_calibration_" + name.replace("(", "").replace(")", "").replace("+", "").replace(" ", "_") + ".png"
    _, res = calibration(test["CRC"].values, test[f"p_{name}"].values, name, fname)
    rows.append({
        "Model": name,
        "AUROC": roc_auc_score(test["CRC"], test[f"p_{name}"]),
        "Cal. intercept": res.loc["Calibration intercept", "Estimate"],
        "Cal. intercept SE": res.loc["Calibration intercept", "Std. Error"],
        "Cal. slope": res.loc["Calibration slope", "Estimate"],
        "Cal. slope SE": res.loc["Calibration slope", "Std. Error"],
        "AIC (train)": models[name].aic,
    })
comparison = pd.DataFrame(rows).set_index("Model")
out("\nComparison on the test set (same split as B-E):")
out(comparison.round(4).to_string())

roc_plot([(name, test["CRC"], test[f"p_{name}"]) for name in formulas],
         "H_roc_comparison.png", "ROC curves of the three models")

# Robustness: repeat over 20 random 70/30 splits
out("\nRobustness check: 20 random 70/30 train-test splits")
results = {name: [] for name in formulas}
for s in range(20):
    r = np.random.default_rng(s)
    idx = r.choice(len(data), size=int(0.7 * len(data)), replace=False)
    tr = data.iloc[idx]
    te = data.drop(index=data.index[idx])
    for name, f in formulas.items():
        p = fit_logit(f, tr).predict(te).values
        pc = np.clip(p, 1e-8, 1 - 1e-8)
        lp = np.log(pc / (1 - pc))
        cm = sm.GLM(te["CRC"].values, sm.add_constant(lp), family=sm.families.Binomial()).fit()
        results[name].append((roc_auc_score(te["CRC"], p), cm.params[0], cm.params[1]))

base_auc = np.array([r[0] for r in results["Baseline"]])
robust = pd.DataFrame({
    name: {
        "AUROC mean": np.mean([r[0] for r in res]),
        "AUROC sd": np.std([r[0] for r in res], ddof=1),
        "Cal. intercept mean": np.mean([r[1] for r in res]),
        "Cal. slope mean": np.mean([r[2] for r in res]),
        "AUROC > baseline (of 20)": int(np.sum(np.array([r[0] for r in res]) > base_auc)),
    }
    for name, res in results.items()
}).T
out(robust.round(4).to_string())

REPORT.close()
print(f"\nAll output saved in: {OUT}")

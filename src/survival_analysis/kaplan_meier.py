import pandas as pd
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter

""" Required packages: numpy, pandas, lifelines
About lifelines KM model: https://lifelines.readthedocs.io/en/latest/Quickstart.html """

# Paths with the splitted_samples files vvv
PATH_STAGE_LOW  = "G:\\CulebraExt\\splitted_samples\\stage_low.csv"
PATH_STAGE_HIGH = "G:\\CulebraExt\\splitted_samples\\stage_high.csv"

# kind of survival data to be considered vvv
SURV = "overall_surv"
# SURV = "dis_free_surv"
# SURV = "dis_specific_surv"

################################################### dataframe preparation
df_low = pd.read_csv(PATH_STAGE_LOW, index_col = 0) # load the dataframe
df_low = df_low[[f"{SURV}_months", SURV]]           # consider only the time and event columns (respectively)
df_low = df_low.dropna()                            # remove rows with missing values

df_high = pd.read_csv(PATH_STAGE_HIGH, index_col = 0)
df_high = df_high[[f"{SURV}_months", SURV]]
df_high = df_high.dropna()

################################################### elaboration of KM curves
kmf = KaplanMeierFitter()

T_low = pd.to_numeric(df_low[f"{SURV}_months"].str.replace(',', '.'))
E_low = df_low[SURV].astype(bool)
kmf.fit(T_low, event_observed = E_low, label = f"Stage Low ({SURV})")
ax = kmf.plot_survival_function()

T_high = pd.to_numeric(df_high[f"{SURV}_months"].str.replace(',', '.'))
E_high = df_high[SURV].astype(bool)
kmf.fit(T_high, event_observed = E_high, label = f"Stage High ({SURV})")
ax = kmf.plot_survival_function(ax = ax)

plt.show()

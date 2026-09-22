# trigger_eff

Trigger-efficiency dataframes, based on the HNL `nopreselect` + `savePfp` pattern:
all PFPs kept, no cosmic / nu_score / FV preselection, GENIE `mcnu` truth matched via `tmatch`.

## Example commands

```bash
# MC (GENIE) — no preselection, save all PFPs
python run_df_maker.py \
  -c ./analysis_village/trigger_eff/configs/trigger_eff_mcnu_nopreselect_savepfp.py \
  -l input.list \
  -o trigger_eff_mc.df

# Data — no preselection, save all PFPs
python run_df_maker.py \
  -c ./analysis_village/trigger_eff/configs/trigger_eff_data_nopreselect_savepfp.py \
  -l input.list \
  -o trigger_eff_data.df
```

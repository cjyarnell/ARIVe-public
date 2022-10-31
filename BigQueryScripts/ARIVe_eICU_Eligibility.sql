-- generates a table with a row for every icu stay_id and columns for each of the eligibility criteria
-- (ie fiO2 40% or more, not yet invasively ventilated, 
-- results saved as ARIVe_eICU_Eligibility

with o2tv as (
select 
patientunitstayid,
charttime,
last_value(fio2 ignore nulls) 
  over (partition by patientunitstayid order by charttime
        range between unbounded preceding and current row)
  as last_fio2,
last_value(o2_device ignore nulls) 
  over (partition by patientunitstayid order by charttime
        range between unbounded preceding and current row)
  as last_o2_device,
case when max(o2_device = 6) over
        (partition by patientunitstayid order by charttime
          range between unbounded preceding and current row)
        then 1 else 0 end as ett_ever,
case when max(o2_device = 7) over
        (partition by patientunitstayid order by charttime
          range between unbounded preceding and 10080 following)
        then 1 else 0 end as trach_ever
from(
select 
patientunitstayid,
charttime,
variable,
value
from alpine-scholar-292916.ARIVe.eICU_oxygenFlowFiO2 o2ff
where (value >= 21 and value <= 100)
union all
select 
patientunitstayid,
charttime,
variable,
value
from alpine-scholar-292916.ARIVe.eICU_oxygenDevices o2dev
) pivot(max(value) for variable IN ("fio2","o2_device"))
order by 1,2
)

select distinct
patientunitstayid,
first_value(charttime) over 
  (partition by patientunitstayid order by charttime) 
  as eligibletime,
first_value(last_fio2) over
  (partition by patientunitstayid order by charttime) 
  as eligiblefio2,
first_value(last_o2_device ignore nulls) over
  (partition by patientunitstayid order by charttime) 
  as baseline_o2_device,
first_value(ett_ever) over
  (partition by patientunitstayid order by charttime) 
  as ett_ever
from (select * from o2tv where last_o2_device is not null)
where --last_fio2 >= 40
   --   last_fio2 > 21
   last_o2_device is not null
  and last_o2_device >= 1
  and ett_ever = 0
  and trach_ever = 0
  and charttime <= 1440


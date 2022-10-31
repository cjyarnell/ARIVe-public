-- ARIVe eICU timevarying

with nc as (
  select 
  a.patientunitstayid,
  nursingchartoffset,
  nursingchartcelltypecat,
  nursingchartcelltypevallabel,
  nursingchartcelltypevalname,
  nursingchartvalue
from alpine-scholar-292916.ARIVe.ARIVe_eICU_Eligibility a
left join physionet-data.eicu_crd.nursecharting b
on a.patientunitstayid = b.patientunitstayid
where nursingchartcelltypecat IN ("Vital Signs","Scores")
), elig as (
  select * from alpine-scholar-292916.ARIVe.ARIVe_eICU_Eligibility
)

select * from(
select
nc.patientunitstayid,
nursingchartoffset-eligibletime as charttime,
case  when nursingchartcelltypevalname = "Respiratory Rate" then "resp_rate"
      when nursingchartcelltypevalname = "O2 Saturation" then "spo2"
      when nursingchartcelltypevalname = "Heart Rate" then "heart_rate"
      when nursingchartcelltypevalname in ("Non-Invasive BP Systolic",
                                           "Invasive BP Systolic") then "sbp"
                                           end as variable,
case when nursingchartcelltypevalname = "GCS Total" then cast(regexp_extract(nursingchartvalue, "[0-9]*") as numeric)
  else cast(nursingchartvalue as numeric) end as value
from nc 
left join elig
on nc.patientunitstayid = elig.patientunitstayid 
where nursingchartcelltypevalname in ("Respiratory Rate",
                                      "O2 Saturation",
                                      "Heart Rate",
                                      "Non-Invasive BP Systolic",
                                           "Invasive BP Systolic")

union all

select
o2.patientunitstayid,
o2.charttime - eligibletime as charttime ,
variable,
value
from elig 
left join alpine-scholar-292916.ARIVe.eICU_oxygenDevices o2
on o2.patientunitstayid = elig.patientunitstayid

union all

select
o2.patientunitstayid,
o2.charttime - eligibletime as charttime,
case when fio2_method = "recorded" then "fio2_recorded"
     when fio2_method = "calculated" then "fio2_calculated" 
     end as variable,
value
from elig 
left join alpine-scholar-292916.ARIVe.eICU_oxygenFlowFiO2 o2  
on o2.patientunitstayid = elig.patientunitstayid

union all

select
elig.patientunitstayid,
chartoffset - eligibletime as charttime,
"pressor" as variable,
greatest(norepinephrine, epinephrine, vasopressin, phenylephrine) as value
from elig
left join physionet-data.eicu_crd_derived.pivoted_infusion prs
on elig.patientunitstayid = prs.patientunitstayid

union all

select
elig.patientunitstayid,
chartoffset - eligibletime as charttime,
"gcs" as variable,
gcs as value
from elig
left join physionet-data.eicu_crd_derived.pivoted_gcs g
on elig.patientunitstayid = g.patientunitstayid

) where charttime <= 24*60*28
order by 1,2,3
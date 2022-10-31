-- ARIVe eICU baseline table

with elig as (
  select * from alpine-scholar-292916.ARIVe.ARIVe_eICU_Eligibility
)

, hospcoeff as (
    SELECT  
  hospitalid,
  avg(case when ethnicity in ("African American","Hispanic","Asian", "Native American") then 1 else 0 end) as prop_nonwhite
  FROM `physionet-data.eicu_crd.patient`
  where ethnicity in ("Caucasian","African American","Hispanic","Asian", "Native American")
  group by 1
  order by 2
)

select distinct
elig.patientunitstayid,
pt.patienthealthsystemstayid,
pt.uniquepid,

pt.gender,
pt.age,
pt.ethnicity,

pt.hospitalid,
pt.wardid,
pt.unittype,

pt.apacheadmissiondx,
pt.admissionheight,
pt.admissionweight,

pt.hospitaladmitoffset,
pt.hospitaladmitsource,

elig.eligibletime,
pt.unitdischargeoffset,
pt.unitdischargelocation,
pt.hospitaldischargeoffset,
pt.hospitaldischargestatus,
pt.hospitaldischargelocation,

hosp.numbedscategory,
hosp.teachingstatus,
hosp.region,

max(case when regexp_contains(lower(pasthistoryvalue), "chf") then 1 else 0 end) as chf,
max(case when regexp_contains(lower(pasthistoryvalue), "copd") then 1 else 0 end) as copd,
max(case when regexp_contains(lower(pasthistoryvalue), "dementia") then 1 else 0 end) as dementia,
max(case when regexp_contains(lower(pasthistorypath), "cancer") then 1 else 0 end) as cancer,

hospcoeff.prop_nonwhite,
pt.unitvisitnumber

from elig
left join physionet-data.eicu_crd.patient pt
on elig.patientunitstayid = pt.patientunitstayid
left join physionet-data.eicu_crd.hospital hosp
on hosp.hospitalid = pt.hospitalid
left join physionet-data.eicu_crd.pasthistory pmh
on elig.patientunitstayid = pmh.patientunitstayid
left join hospcoeff 
on pt.hospitalid = hospcoeff.hospitalid
--where pasthistoryoffset <= eligibletime
group by 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,28, 29
order by 1

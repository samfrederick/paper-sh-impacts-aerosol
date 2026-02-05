#!/bin/bash

home_path=$(pwd)
processed_data_dir="../data"
processed_path="${home_path}/${processed_data_dir}"
#          no het   low het  med het   high het
slurm_ids=(224670   224841   2726708   231241)

#----------------------------------------------------------------------------
# DATA FOR FIGURES 3,4,6,7,10,12

time_idx=36
#z_idx=45
aero_variables="nh3,hno3,oh,pmc_NH4,pmc_NO3,pmc_SO4,ccn_001,ccn_003,ccn_006,ccn_010"
wrf_variables="ALT"

scenario_number=0
for id in "${slurm_ids[@]}"
do 
    scenario="scenario-${scenario_number}"
    if [[ "$scenario" == "scenario-0" ]]; then
        scenario="no-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-1" ]]; then
        scenario="low-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-2" ]]; then
        scenario="medium-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-3" ]]; then
        scenario="high-heterogeneity"
    fi

    if [[ "$scenario" == "no-heterogeneity" ]]; then
        path="/data/keeling/a/sf20/d/les_output/wrf-partmc/slurm-${id}"
    else
        path="/data/keeling/a/sf20/e/wrf-partmc-output/slurm-${id}"
    fi

    cd ${path}
    aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
    subset_filename=${scenario}_subset_t${time_idx}.nc

    output_path="${processed_path}/${subset_filename}"

    if [ ! -f "$output_path" ]; then
        echo -e "Creating subset of ${scenario}"
        ncks -v ${aero_variables} -d Time,${time_idx} ${aero_data} ${output_path}
        wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
        ncks -A -v ${wrf_variables} -d Time,${time_idx} ${wrf_data} ${output_path}
        size=$(du -sh ${processed_path}/${subset_filename})
        echo -e "File size: ${size}\n"
    fi 

    ((scenario_number+=1))
done 

# no het, no ammonia
scenario="no-heterogeneity-no-nh4"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-230894
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
output_path="${processed_path}/${scenario}_subset_t${time_idx}.nc"
if [ ! -f "$output_path" ]; then
    echo -e "Creating subset of ${scenario}"
    ncks -v ${aero_variables} -d Time,${time_idx} ${aero_data} ${output_path}
    wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
    ncks -A -v ${wrf_variables} -d Time,${time_idx} ${wrf_data} ${output_path}
    size=$(du -sh ${processed_path}/${scenario}_subset_t${time_idx}.nc)
    echo -e "File size: ${size}\n" 
fi

# high het, no ammonia
scenario="high-heterogeneity-no-nh4"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-2727876
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
output_path="${processed_path}/${scenario}_subset_t${time_idx}.nc"
if [ ! -f "$output_path" ]; then
    echo -e "Creating subset of ${scenario}"
    ncks -v ${aero_variables} -d Time,${time_idx} ${aero_data} ${output_path}
    wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
    ncks -A -v ${wrf_variables} -d Time,${time_idx} ${wrf_data} ${output_path}
    size=$(du -sh ${processed_path}/${scenario}_subset_t${time_idx}.nc)
    echo -e "File size: ${size}\n" 
fi


#----------------------------------------------------------------------------
# DATA FOR FIGURE 5
time_idx=36
z_idx=60

num_vars=$(printf "num_a%03d," {1..100})
num_vars=${num_vars%,}
mass_vars=$(printf "mass_a%03d," {1..100})
mass_vars=${mass_vars%,}
aerodist_variables="${num_vars},${mass_vars}"
aerobin_variables="BIN_EDGES,BIN_CENTERS"

scenario_number=0
for id in "${slurm_ids[@]}"
do 
    scenario="scenario-${scenario_number}"
    if [[ "$scenario" == "scenario-0" ]]; then
        scenario="no-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-1" ]]; then
        scenario="low-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-2" ]]; then
        scenario="medium-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-3" ]]; then
        scenario="high-heterogeneity"
    fi

    if [[ "$scenario" == "no-heterogeneity" ]]; then
        path="/data/keeling/a/sf20/d/les_output/wrf-partmc/slurm-${id}"
    else
        path="/data/keeling/a/sf20/e/wrf-partmc-output/slurm-${id}"
    fi

    cd ${path}
    aerodist_data=${path}/aerosol_dist_d01_2023-03-20_09:00:00
    subset_filename=${scenario}_size-dist_subset_t0_t${time_idx}_z${z_idx}.nc

    output_path="${processed_path}/${subset_filename}"
    if [ ! -f "$output_path" ]; then
        echo -e "Creating size distribution subset of ${scenario}"
        ncks -v ${aerodist_variables} -d Time,0 -d Time,${time_idx} -d bottom_top,${z_idx} ${aerodist_data} ${output_path}
        ncks -A -v ${aerobin_variables} -d Time,0 -d Time,${time_idx} ${aerodist_data} ${output_path}
        size=$(du -sh ${processed_path}/${subset_filename})
        echo -e "File size: ${size}\n"
    fi

    ((scenario_number+=1))
done 

#----------------------------------------------------------------------------
# Cross section averaged quantities (all t, all z)
# Figure 11, etc?

aero_variables="ccn_001,ccn_003,ccn_006,ccn_010,pmc_NH4,pmc_NO3,pmc_SO4,nh3,hno3"
avg_dimensions="south_north,west_east" # average for each domain horizontal cross section

scenario_number=0
for id in "${slurm_ids[@]}"
do 
    scenario="scenario-${scenario_number}"
    if [[ "$scenario" == "scenario-0" ]]; then
        scenario="no-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-1" ]]; then
        scenario="low-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-2" ]]; then
        scenario="medium-heterogeneity"
    fi
    if [[ "$scenario" == "scenario-3" ]]; then
        scenario="high-heterogeneity"
    fi

    if [[ "$scenario" == "no-heterogeneity" ]]; then
        path="/data/keeling/a/sf20/d/les_output/wrf-partmc/slurm-${id}"
    else
        path="/data/keeling/a/sf20/e/wrf-partmc-output/slurm-${id}"
    fi

    cd ${path}
    aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
    subset_filename=${scenario}_time-height-avg_subset.nc

    output_path="${processed_path}/${subset_filename}"
    if [ ! -f "$output_path" ]; then
        echo -e "Creating time-height subset of ${scenario}"
        ncks -v ${aero_variables} ${aero_data} ${output_path}
        ncwa -O -a ${avg_dimensions} -v ${aero_variables} ${output_path} ${output_path}

        size=$(du -sh ${processed_path}/${subset_filename})
        echo -e "File size: ${size}\n"
    fi

    ((scenario_number+=1))
done 

# no het no emiss. after t = 4 h
scenario="no-heterogeneity-emis-off-t4"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-229754
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
subset_filename=${scenario}_time-height-avg_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating time-height subset of ${scenario}"
    ncks -v ${aero_variables} ${aero_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${aero_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n"
fi

# high het no emiss. after t = 4 h
scenario="high-heterogeneity-emis-off-t4"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-230150
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
subset_filename=${scenario}_time-height-avg_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating time-height subset of ${scenario}"
    ncks -v ${aero_variables} ${aero_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${aero_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n"
fi 

# No het, high RH
scenario="no-heterogeneity-high-RH"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-234773
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
subset_filename=${scenario}_time-height-avg_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating time-height subset of ${scenario}"
    ncks -v ${aero_variables} ${aero_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${aero_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n" 
fi 

# High het, high RH
scenario="high-heterogeneity-high-RH"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-234798
cd ${path}
aero_data=${path}/aerosols_d01_2023-03-20_09:00:00
subset_filename=${scenario}_time-height-avg_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating time-height subset of ${scenario}"
    ncks -v ${aero_variables} ${aero_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${aero_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n" 
fi 

# Data for Figure 2 (Temp, RH time height plots)
wrf_variables="QVAPOR,P,PB,T"
avg_dimensions="south_north,west_east" # average for each domain horizontal cross section
scenario="no-heterogeneity"
path="/data/keeling/a/sf20/d/les_output/wrf-partmc/slurm-${slurm_ids[0]}"

cd ${path}
wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
subset_filename=${scenario}_met-vars_subset.nc

output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating met subset of ${scenario}"
    ncks -v ${wrf_variables} ${wrf_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${wrf_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n"
fi 


wrf_variables="QVAPOR,P,PB,T"
avg_dimensions="south_north,west_east" # average for each domain horizontal cross section
scenario="no-heterogeneity-high-RH"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-234773
cd ${path}
wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
subset_filename=${scenario}_met-vars_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating met subset of ${scenario}"
    ncks -v ${wrf_variables} ${wrf_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${wrf_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n"
fi 

# TODO add met data for high het high rh scenario
wrf_variables="QVAPOR,P,PB,T"
avg_dimensions="south_north,west_east" # average for each domain horizontal cross section
scenario="high-heterogeneity-high-RH"
path=/data/keeling/a/sf20/e/wrf-partmc-output/slurm-234798
cd ${path}
wrf_data=${path}/wrfout_d01_2023-03-20_09:00:00
subset_filename=${scenario}_met-vars_subset.nc
output_path="${processed_path}/${subset_filename}"
if [ ! -f "$output_path" ]; then
    echo -e "Creating met subset of ${scenario}"
    ncks -v ${wrf_variables} ${wrf_data} ${output_path}
    ncwa -O -a ${avg_dimensions} -v ${wrf_variables} ${output_path} ${output_path}
    size=$(du -sh ${processed_path}/${subset_filename})
    echo -e "File size: ${size}\n"
fi 
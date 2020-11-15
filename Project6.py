# Homework #6 

# 1. Find a way to represent the entire time series in a single file. 
# You may either save a visualization for every data on separate pages in a pdf or you may create an MP4 or GIF file. 
# Make sure that for every frame, the correct date is indicated. Differentiate your answer from the template provided in at least 3 ways. 
# Identify how your output is differentiated.

# 2. Select a different time series data (e.g., unemployment, incomes, etc...) set at the county level that has at least 20 observations. 
# Create a map of the latest obervation of the data. Like with the covid data in (1), create a file that holds a visualization for every observation.

# ECON 611 Only

# 3. Download and plot data relating to the magnitude of the policy response by U.S. state or by all nations across the globe. 
# Plot the data on a map as it changes over time using an MP4 file or a GIF.


#createCOVID19StateAndCountyVisualization.py
import geopandas
import numpy as np
import pandas as pd
# We won't actually use datetime directly. Since the dataframe index will use 
# data formatted as datetime64, I import it in case I need to use the datetime
# module to troubleshoot later 
import datetime
# you could technically call many of the submodules from matplotlib using mpl., 
#but for convenience we explicitly import submodules. These will be used for 
# constructing visualizations
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mtick
import datadotworld as dw


def import_geo_data(filename, index_col = "Date", FIPS_name = "FIPS"):
    # import county level shapefile
    map_data = geopandas.read_file(filename = filename,                                   
                                   index_col = index_col)
    # rename fips code to match variable name in COVID-19 data
    map_data.rename(columns={"State":"state"},
                    inplace = True)
    # Combine statefips and county fips to create a single fips value
    # that identifies each particular county without referencing the 
    # state separately
    map_data[FIPS_name] = map_data["STATEFP"].astype(str) + \
        map_data["COUNTYFP"].astype(str)
    map_data[FIPS_name] = map_data[FIPS_name].astype(np.int64)
    # set FIPS as index
    map_data.set_index(FIPS_name, inplace=True)
    
    return map_data

def import_covid_data(filename, FIPS_name):
    # Load COVID19 county data using datadotworld API
    # Data provided by Johns Hopkins, file provided by Associated Press
    dataset = dw.load_dataset("associatedpress/johns-hopkins-coronavirus-case-tracker",
                              auto_update = True)
    # the dataset includes multiple dataframes. We will only use #2
    covid_data = dataset.dataframes["2_cases_and_deaths_by_county_timeseries"]
    # Include only oberservation for political entities within states
    # i.e., not territories, etc... drop any nan fip values with covid_data[FIPS_name] > 0
    covid_data = covid_data[covid_data[FIPS_name] < 57000]
    covid_data = covid_data[covid_data[FIPS_name] > 0]

    # Transform FIPS codes into integers (not floats)
    covid_data[FIPS_name] = covid_data[FIPS_name].astype(int)
    covid_data.set_index([FIPS_name, "date"], inplace = True)
    # Prepare a column for state abbreviations. We will draw these from a
    # dictionary created in the next step.
    covid_data["state_abr"] = ""
    for state, abr in state_dict.items():
        covid_data.loc[covid_data["state"] == state, "state_abr"] = abr
    # Create "Location" which concatenates county name and state abbreviation 
    covid_data["Location"] = covid_data["location_name"] + ", " + \
        covid_data["state_abr"]

    return covid_data

def create_covid_geo_dataframe(covid_data, map_data, dates):
    # create geopandas dataframe with multiindex for date
    # original geopandas dataframe had no dates, so copies of the df are 
    # stacked vertically, with a new copy for each date in the covid_data index
    #(dates is a global)
    i = 0
    for date in dates:
        # select county observations from each date in dates
        df = covid_data[covid_data.index.get_level_values("date")==date]
        # use the fips_codes from the slice of covid_data to select counties
        # from the map_data index,making sure that the map_data index matches
        # the covid_data index
        counties = df.index.get_level_values("fips_code")
        agg_df = map_data.loc[counties]
        # each row for agg_df will reflect that 
        agg_df["date"] = date
        if i == 0:
            # create the geodataframe, select coordinate system (.crs) to
            # match map_data.crs
            matching_gpd = geopandas.GeoDataFrame(agg_df, crs = map_data.crs)
            i += 1
        else:
            # after initial geodataframe is created, stack a dataframe for
            # each date in dates. Once completed, index of matching_gpd
            # will match index of covid_data
            matching_gpd = matching_gpd.append(agg_df, ignore_index = False)         
    # Set mathcing_gpd index as["fips_code", "date"], liked covid_data index
    matching_gpd.reset_index(inplace=True)
    matching_gpd.set_index(["fips_code","date"], inplace = True)
    # add each column from covid_data to mathcing_gpd
    for key, val in covid_data.items():
        matching_gpd[key] = val

    return matching_gpd       

def create_new_vars(covid_data, moving_average_days):
    # use a for loop that performs the same operations on data for cases and for deaths
    for key in ["cases", "deaths"]:
        # create a version of the key with the first letter capitalized
        cap_key = key.title()
        covid_data[cap_key + " per Million"] = covid_data["cumulative_" + key]\
            .div(covid_data["total_population"]).mul(10 ** 6)
        # generate daily data normalized per million population by taking the daily difference within each
        # entity (covid_data.index.names[0]), dividing this value by population and multiplying that value
        # by 1 million 10 ** 6
        covid_data["Daily " + cap_key + " per Million"] = \
            covid_data["cumulative_" + key ].groupby(covid_data.index.names[0])\
            .diff(1).div(covid_data["total_population"]).mul(10 ** 6)
        # taking the rolling average; choice of number of days is passed as moving_average_days
        covid_data["Daily " + cap_key + " per Million MA"] = covid_data["Daily " + \
                  cap_key + " per Million"].rolling(moving_average_days).mean()

# I include this dictionary to convenienlty cross reference state names and
# state abbreviations.
state_dict = {
    'Alabama': 'AL', 'Alaska': 'AK', 'Arizona': 'AZ',
    'Arkansas': 'AR', 'California': 'CA', 'Colorado': 'CO', 'Connecticut': 'CT', 
    'Delaware': 'DE', 'District of Columbia': 'DC', 'Florida': 'FL', 
    'Georgia': 'GA', 'Hawaii': 'HI', 'Idaho': 'ID', 'Illinois': 'IL',
    'Indiana': 'IN', 'Iowa': 'IA','Kansas': 'KS', 'Kentucky': 'KY',
    'Louisiana': 'LA', 'Maine': 'ME', 'Maryland': 'MD', 'Massachusetts': 'MA',
    'Michigan': 'MI', 'Minnesota': 'MN', 'Mississippi': 'MS', 'Missouri': 'MO',
    'Montana': 'MT', 'Nebraska': 'NE', 'Nevada': 'NV', 'New Hampshire': 'NH',
    'New Jersey': 'NJ', 'New Mexico': 'NM', 'New York': 'NY', 'North Carolina': 'NC',
    'North Dakota': 'ND', 'Ohio': 'OH', 'Oklahoma': 'OK',
    'Oregon': 'OR', 'Pennsylvania': 'PA', 'Rhode Island': 'RI',
    'South Carolina': 'SC', 'South Dakota': 'SD', 'Tennessee': 'TN', 'Texas': 'TX',
    'Utah': 'UT', 'Vermont': 'VT', 'Virginia': 'VA',
    'Washington': 'WA', 'West Virginia': 'WV', 'Wisconsin': 'WI', 'Wyoming': 'WY'}

plt.rcParams['axes.ymargin'] = 0
plt.rcParams['axes.xmargin'] = 0
plt.rcParams.update({'font.size': 32})

if "data_processed" not in locals():
    fips_name = "fips_code"
    covid_filename = "COVID19DataAP.csv"
    # rename_FIPS matches map_data FIPS with COVID19 FIPS name
    map_data = import_geo_data(filename = "countiesWithStatesAndPopulation.shp",
                    index_col = "Date", FIPS_name= fips_name)
    covid_data = import_covid_data(filename = covid_filename, FIPS_name = fips_name)
    # dates will be used to create a geopandas DataFrame with multiindex 
    dates = sorted(list(set(covid_data.index.get_level_values("date"))))
    covid_data = create_covid_geo_dataframe(covid_data, map_data, dates)
    moving_average_days = 7
    create_new_vars(covid_data, moving_average_days)
    start_date = "03-15-2020"     
    end_date = dates[-1]
    # once data is processed, it is saved in the memory
    # the if statement at the top of this block of code instructs the computer
    # not to repeat these operations 
    data_processed = True
    

def select_data_within_bounds(data, minx, miny, maxx, maxy):
    data = data[data.bounds["maxx"] <= maxx]
    data = data[data.bounds["maxy"] <= maxy]
    data = data[data.bounds["minx"] >= minx]
    data = data[data.bounds["miny"] >= miny]
    
    return data

date = dates[-1]

if "map_bounded" not in locals():
    minx = -127
    miny = 23
    maxx = -58
    maxy = 54
    covid_map_data = select_data_within_bounds(covid_data, minx, miny, maxx, maxy)
    map_bounded = True


fig, ax = plt.subplots(figsize=(18,8),
        subplot_kw = {'aspect': 'equal'})   
plt.rcParams.update({"font.size": 30})
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
key = "Deaths per Million"
# this time we replace 0 values with 1
# so that these values show up as beige instead of as white
# when color axis is logged
df = covid_map_data[covid_map_data.index.get_level_values("date")==date].replace(0,1)
# set range of colorbar
vmin = 1 
vmax = df[key].max()
# choose colormap
cmap = cm.get_cmap('Blue', 4)
# format colormap
norm = cm.colors.LogNorm(vmin = vmin, vmax = vmax)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# empty array for the data range
sm._A = []
# prepare space for colorbar
divider = make_axes_locatable(ax)
size = "10%" 
cax = divider.append_axes("right", size = size, pad = 0.1)
# add colorbar to figure
cbar = fig.colorbar(sm, cax=cax, cmap = cmap)
cbar.ax.tick_params(labelsize=20)
vals = list(cbar.ax.get_yticks())
vals.append(vmax)
# format colorbar values as int
cbar.ax.set_yticklabels([int(x) for x in vals])
cbar.ax.set_ylabel(key, fontsize = 20)


df.plot(ax=ax, cax = cax, column=key, vmin=vmin ,vmax = vmax, 
             cmap = cmap, legend=False, linewidth=.75, edgecolor='lightgrey', 
             norm = norm)
ax.set_title(str(date)[:10] + "\n" + "COVID-19 in the U.S.", fontsize = 30)

pp = PdfPages('Project#6.pdf')
pp.savefig(map_data)
pp.savefig(covid_data)
pp.close()

#Change formatting in 3 different ways:
    # I changed the map color to blue
    # I made the color bar bigger in order to see it more clear
    # I increased line width to make it easier to see
    

#Problem #2


def import_geo_data(filename, index_col = "Date", FIPS_name = "FIPS"):
    # import county level shapefile
    map_data = geopandas.read_file(filename = filename,                                   
                                   index_col = index_col)
    # rename fips code to match variable name in COVID-19 data
    map_data.rename(columns={"State":"state"},
                    inplace = True)
    # Combine statefips and county fips to create a single fips value
    # that identifies each particular county without referencing the 
    # state separately
    map_data[FIPS_name] = map_data["STATEFP"].astype(str) + \
        map_data["COUNTYFP"].astype(str)
    map_data[FIPS_name] = map_data[FIPS_name].astype(np.int64)
    # set FIPS as index
    map_data.set_index(FIPS_name, inplace=True)
    
    return map_data


def import_Umemployment_data(filename, FIPS_name):
    # data provided by USDA
    # https://www.ers.usda.gov/data-products/county-level-data-sets/download-data/
    dataset = dw.load_dataset("unemployment.csv",
                              auto_update = True)
    # the dataset includes multiple dataframes. We will only use
    unemployment_data = dataset.dataframes["Unemployment_rate_2019"]
    # drop any nan fip values with covid_data[FIPS_name] > 0
    unemployment_data = unemployment_data[covid_data[FIPS_name] < 57000]
    unemployment_data = unemployment_data[covid_data[FIPS_name] > 0]

    # Transform FIPS codes into integers (not floats)
    unemployment_data[FIPS_name] = unemployment_data[FIPS_name].astype(int)
    unemployment_data.set_index([FIPS_name, "date"], inplace = True)
    # Prepare a column for state abbreviations. We will draw these from a
    # dictionary created in the next step.
    unemployment_data["state_abr"] = ""
    for state, abr in state_dict.items():
        unemployment_data.loc[unemployment_data["state"] == state, "state_abr"] = abr
    # Create "Location" which concatenates county name and state abbreviation 
    unemployment_data["Location"] = unemployment_data["location_name"] + ", " + \
        unemployment_data["state_abr"]

    return unemployment_data

def create_unemployment_geo_dataframe(unemployment_data, map_data, dates):
    # create geopandas dataframe with multiindex for date
    # original geopandas dataframe had no dates, so copies of the df are 
    # stacked vertically, with a new copy for each date in the covid_data index
    #(dates is a global)
    i = 0
    for date in dates:
        # select county observations from each date in dates
        df = unemployment_data[unemployment_data.index.get_level_values("date")==date]
        # use the fips_codes from the slice of covid_data to select counties
        # from the map_data index,making sure that the map_data index matches
        # the covid_data index
        counties = df.index.get_level_values("fips_code")
        agg_df = map_data.loc[counties]
        # each row for agg_df will reflect that 
        agg_df["date"] = date
        if i == 0:
            # create the geodataframe, select coordinate system (.crs) to
            # match map_data.crs
            matching_gpd = geopandas.GeoDataFrame(agg_df, crs = map_data.crs)
            i += 1
        else:
            # after initial geodataframe is created, stack a dataframe for
            # each date in dates. Once completed, index of matching_gpd
            # will match index of covid_data
            matching_gpd = matching_gpd.append(agg_df, ignore_index = False)         
    # Set mathcing_gpd index as["fips_code", "date"], liked covid_data index
    matching_gpd.reset_index(inplace=True)
    matching_gpd.set_index(["fips_code","date"], inplace = True)
    # add each column from covid_data to mathcing_gpd
    for key, val in covid_data.items():
        matching_gpd[key] = val

    return matching_gpd 

pp1 = PdfPages('Project#6.pdf')
pp1.savefig(unemployment_data)
pp1.close()


#Problem #3


df_Policy = datareader(OxCGRT_latest.csv)
# set range of colorbar
vmin = 1 
vmax = df[key].max()
# choose colormap
cmap = cm.get_cmap('Blue', 4)
# format colormap
norm = cm.colors.LogNorm(vmin = vmin, vmax = vmax)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# empty array for the data range
sm._A = []
# prepare space for colorbar
divider = make_axes_locatable(ax)
size = "5%" 
cax = divider.append_axes("right", size = size, pad = 0.1)
# add colorbar to figure
cbar = fig.colorbar(sm, cax=cax, cmap = cmap)
cbar.ax.tick_params(labelsize=18)
vals = list(cbar.ax.get_yticks())
vals.append(vmax)
# format colorbar values as int
cbar.ax.set_yticklabels([int(x) for x in vals])
cbar.ax.set_ylabel(key, fontsize = 20)


df.plot(ax=ax, cax = cax, column=key, vmin=vmin ,vmax = vmax, 
             cmap = cmap, legend=False, linewidth=.5, edgecolor='lightgrey', 
             norm = norm)
ax.set_title(str(date)[:10] + "\n" + "C1_flag", fontsize = 30)


if "map_bounded" not in locals():
    minx = -127
    miny = 23
    maxx = -58
    maxy = 54
    covid_data = select_data_within_bounds(covid_data, minx, miny, maxx, maxy)
    map_bounded = True


for key in keys:
    log = False if "Daily" in key else True
    # this time we replace 0 values with 1
    # so that these values show up as beige  instead of as white
    # when color axis is logged
    vmin = 1 if log else 0 
    vmax = df[key][df.index.get_level_values("date") == date].max()
    # Create colorbar as a legend
    cmap = cm.get_cmap('Reds', 4)
    if log:
        norm = cm.colors.LogNorm(vmin=vmin, vmax =vmax)
    else:
        norm = cm.colors.Normalize(vmin = vmin, vmax = vmax)

    fig, ax = plt.subplots(figsize=(18,8),
        subplot_kw = {'aspect': 'equal'})   
    plt.rcParams.update({"font.size": 30})
    plt.xticks(fontsize = 25)
    plt.yticks(fontsize = 25)
    # the functions will unpack the tuple. The same names variable names
    # are used in the function
    kwargs = (df, key, log, date, fig, ax, cmap, norm, vmin, vmax)
    
    plt.show()    

    plt.close()







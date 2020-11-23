import pandas as pd
import pandas_datareader.data as web
import datetime

def gather_data1(stocks, start, end = datetime.datetime.today(), freq = "M"):
    i = 0
    # dct.items() calls key and value that key points to
    for key, val in stocks.items():
        if i == 0:
            # Create dataframe for first variable, then rename column
            df = web.DataReader(val, "yahoo", start, end).resample(freq).mean()
            df.rename(columns = {val:key}, inplace = True) 
            # setting i to None will cause the next block of code to execute,
            # placing data within df instead of creating a new dataframe for
            # each variable
            i = None
        else:
            # If dataframe already exists, add new column
            df[key] = web.DataReader(val, "yahoo", start, end).resample(freq).mean()

    return df
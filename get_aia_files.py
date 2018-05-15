import pandas as pd
from aia_mkmovie import aia_download_files as adf
from multiprocessing import Pool

def wrap_get_files(args):
    print(args)
    return get_files(*args)

def get_files(stime,etime,archive = 'aia_arch/'):
    fobj = adf.download_files(stime,etime,['193','304','335'],'30m',email='jakub.prchlik@cfa.harvard.edu',odir=archive)
    try:
        fobj.get_drms_files()
    except:
        return

def format_dates(data):
    x = data.str.replace('-','/')
    x = x.str.replace('T',' ')
    x = x.str[:-4]

    return x

#read in csv file
aia = pd.read_csv('Sigmoids2007to2017.csv')


#format start times and end times
aia['fmt_stime'] = format_dates(aia.SIG_START)
aia['fmt_etime'] = format_dates(aia.SIG_END)


#remove rows without times
aia.SIG_START.fillna('',inplace=True)
aia.loc[aia.SIG_START != '',:]

#create list for parallel processing
inpt =  aia[['fmt_stime','fmt_etime']].values.tolist()


pool = Pool(processes=4)
out  = pool.map(wrap_get_files,inpt)
pool.close()
pool.join()
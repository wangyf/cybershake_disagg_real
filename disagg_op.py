import argparse
import csv
import subprocess
import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from pathlib2 import Path
import mysql.connector
import sys

import seaborn as sns; 
sns.set(context='paper', style='whitegrid', font_scale=1.3, rc={"xtick.bottom" : True, "ytick.left" : True})
sns.set_style({'legend_frameon':False})


def call_sub(cmd, shell=False, timeout=200):
    p = subprocess.run(cmd, shell=shell, capture_output=True,
                       timeout=timeout, check=True)
    return p


def distance(lon1, lat1, lon2=-117.932587, lat2=33.918633):
    lat1, lon1, lat2, lon2 = np.radians((lat1, lon1, lat2, lon2))
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    d = 0.5 - np.cos(dlat) / 2 + np.cos(lat1) * np.cos(lat2)  * (1 - np.cos(dlon)) / 2
    return 12742 * np.arcsin(np.sqrt(d))


def split_arg(string):
     """Split multiple input args, seperated by ","
     """
     return list(map(float, re.findall(r'\d+\.*\d*', string))) if ',' in string else [int(string)]


ID_map = {
    'Run: ID' : "run_id",
    'Site_ID' : 'site_id',
    'ERF_ID' : 'erf_id',
    'SGT Var ID': 'rup_id',
    'Rup Var Scen ID' : 'rup_var_id',
    'Vel Model ID' : 'vel_model_id',
}
    
def query_ids(line):
    res = {}
    for ID in ('Run: ID', 'Site_ID', 'ERF_ID', 'SGT Var ID',
               'Rup Var Scen ID', 'Vel Model ID'):
        res[ID_map[ID]] = int(re.findall(f'{ID}: (\d+),', line)[0])
    return res


class Disagg:
    def __init__(self, STUDY_ID=7):
        # Start a connection
        self.conn = mysql.connector.connect(
            host="moment.usc.edu",
            user="cybershk_ro",
            password="CyberShake2007",
            database="CyberShake",
            use_unicode=True,
        )
        self.STUDY_ID = STUDY_ID
        self.read_params()
        
        self.query_all_sites()
        self.select_sites_by_dist()
        
        print(f'Location = ({self.lon:.2f}, {self.lat:.2f})')
        print(f'Periods = {self.args["period"]}; Radius = {self.args["radius"]} km')
        print(f'\nOutput to directory:\n\t{self.directory}\n')
        print(f'Neaby sites:\n\t{", ".join(self.sites["name"])}')

        self.run_gmt_plot()

        periods = self.args['period']
        for period in periods:
            dfs, df_sources = self.run(period=period)
            if self.args['plot_sa_map']:
                for i in range(self.args['topk']):
                    self.plot_SA_map(**df_sources.iloc[i].to_dict(), period=period)

    def query_all_sites(self):
        """Query site table
            STUDY_ID = 7 # Study 15.12, with high-frequency statistical implementation
        """
        cmd = (
            f'select S.CS_Short_Name, S.CS_Site_Name, S.CS_Site_Lat, S.CS_Site_Lon '
            f'from CyberShake_Sites S, CyberShake_Runs R '
            f'where R.Study_ID={self.STUDY_ID} '
            f'and R.Site_ID=S.CS_Site_ID;'
        )
        self.all_sites = pd.read_sql(cmd, self.conn)
        self.all_sites.columns = ["name", "deccription", "lat", "lon"]
        self.cursor = self.conn.cursor()

    def read_params(self):
        parser = argparse.ArgumentParser(description="Read input from the command line", usage="For example:\n\tpython disagg.py --lon=-119.2 --lat=34.2 --radius=10.0 period=0.2\n\tpython disagg.py -h\nAt least one of the 'spatial' related parameters are needed, otherwise only GMPE plots are generated.")
        parser.add_argument('--lon', '-lon', default=argparse.SUPPRESS, type=float, help='Longitude at the target location')
        parser.add_argument('--lat', '-lat', default=argparse.SUPPRESS, type=float, help='Latitude at the target location')
        parser.add_argument('--site_name', '-s', type=str, help='Site name of the target location, if lon/lat not provided')
        parser.add_argument('--radius', '-rd', default=10, type=float, help='The radius to search near the target location, default = 10')
        parser.add_argument('--period', '-p', default=[1], type=split_arg, help='Period to disaggregate, default = 1')
        parser.add_argument('--k', '-k', default=1, dest='topk', type=int, help='The number of the most severe sources to be plotted, default = 1')
        parser.add_argument('--spacing', '-sp', default=0.1, type=float, help='Resolution of the SA maps, in degree, default = 0.1')
        parser.add_argument('--spatial_corr_fields', '-scf', default=20, type=int, help='Number of realizations for calculating the spatial correlation')
        parser.add_argument('--spatial_corr_rand', '-scr', type=int, default=argparse.SUPPRESS, help='Random seed for spatially correlated random fields, otherwise unique each time')
        parser.add_argument('--spatial_corr_sigma', '-scs', type=float, default=argparse.SUPPRESS, help='Sigma for spatially correlated random fields, default = 0.6')
        parser.add_argument('--download-interpolated', '-d', type=int, default=argparse.SUPPRESS, help='Whether to download interpoalted grid values, 1 [True] or 0 [False]')
        parser.add_argument('--gmpe', '-g', default='NGAWest_2014_AVG_NOIDRISS', type=str, help='GMPE used for basemap')
        parser.add_argument('--prob', '-b', default=4e-4, type=float, help='Exceedance probability')
        parser.add_argument('--plot_sa_map', '-sa', type=bool, default=True, help='Whether to plot SA maps for random realizations, 1 [True] or 0 [False]')
        parser.add_argument('--study', '-st', default="STUDY_15_12", type=str, help='CyberShake study ID, default = STUDY_15_12')
        parser.add_argument('--cmap', '-c', type=str, default='inferno', help='Colormap of theGMT location map')
        parser.add_argument('--force_update', '-f', default=0, type=bool, help='Whether to force updating queried results, instead of using existing ones directly')
        parser.add_argument('--resolution', '-r', default=15, type=int, help='Resolution of the GMT location map, in arcsec, default=15')
        parser.add_argument('--timeout', '-m', default=2000, type=float, help='Timeout value in case subprocess does not end')
        parser.add_argument('--output', '-o', type=str, help='output directory')
        self.args = vars(parser.parse_args())
        

    def select_sites_by_dist(self, **kwargs):
        """Select sites within {radius} km to (lon, lat)"""
        args = {**self.args, **kwargs}
        sites = self.all_sites
        if self.args['site_name']:
            site = sites[sites['name'] == self.args['site_name']].iloc[0]
            self.lon, self.lat = site.lon, site.lat
        else:
            self.lon, self.lat = args['lon'], args['lat']
        sites['dist'] = sites.apply(
            lambda row: distance(row.lon, row.lat, self.lon, self.lat), axis=1)
        sites.sort_values('dist', inplace=True)
        self.sites = sites[sites['dist'] < args['radius']]
        self.directory = Path(args['output'] or "loc_" + "_".join([str(self.lon), str(self.lat),
            *list(map(str, args["period"]))]))
        self.sites[['lon', 'lat', 'name']].to_csv(Path(self.directory, 'sites_location.txt'),
            sep= ' ', header=None, index=None)
        print(f'\nNearby sites will be written to file: {Path(self.directory, "sites_location.txt")}\n')


    
    def get_disagg(self, file_name):
        with open(file_name, 'r') as fid:
            for i, line in enumerate(fid):
                if i == 13:
                    ids = query_ids(line)

        tmp = pd.read_csv(
            file_name,
            delimiter='\t',
            skiprows=20,
            nrows=1,
        )
        #columns = tmp.columns.str.lower().str.replace(' ', '').to_list() + ['total']
        columns = tmp.columns.str.lower().str.replace(' ', '').to_list()
        df = pd.read_csv(
            file_name,
            delimiter='\t',
            skiprows=22,
            nrows=110,
            names=columns,
        )
        
        columns = ['source_id', 'Contribution', 'ExceedRate', 'source_desc']
        df_sources = pd.read_csv(
            file_name,
            delimiter='\t',
            skiprows=136,
            names=columns,
            usecols=[0, 1, 2, 3],
        )
        
        for ID, v in ids.items():
            df_sources[ID] = v
        return df, df_sources


    def plot_disagg(self, df=None, file_name="", title="", ax=None):
        if df is None:
            df, _ = self.get_disagg(file_name)

        # Specify colors if length known
        colors = mpl.cm.bwr_r(np.linspace(0, 1, len(df.columns) - 3))
        

        # prepare 3d axes
        if ax is None:
            fig = plt.figure(figsize=(10,6), dpi=600)
            ax = Axes3D(fig)

        N = len(df.columns)
        
        # thickness of the bars
        dx, dy = 10, .5

        for i, row in df.iterrows():
            # set up positions for the bars 
            if row.iloc[-1] == 0:
                continue
            xpos = [row.iloc[0]] * (N - 3)
            ypos = [row.iloc[1]] * (N - 3)
            dz = row.iloc[2:-1]
            zpos = dz.cumsum().to_list()
            zpos = [0] + zpos[:-1]
            ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors)


        # label the axes
        ax.set(
            xlabel = 'Rupture Distance (km)',
            ylabel = 'Magnitude',
            zlabel = '%Contribution',
            xlim = (df['dist'].min(), df['dist'].max()),
            ylim = (df['mag'].min(), df['mag'].max()),
        )
        if title:
            ax.set_title(title)
        ax.xaxis.set_major_locator(MultipleLocator(15.))
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.zaxis._axinfo['juggled'] = (1, 2, 0)
        ax.tick_params(axis='both', pad=0)
        #plt.show()


    def run_gmt_plot(self,**kwargs):
        """
        Generate a map of nearby stations at the target location

        Output
        ------
            loc_{lon}_{lat}_{period}/sites_location.pdf

        Example
        -------
        >> Disagg().run_gmt_plot(site_name='LADT')
        >> Disagg().run_gmt_plot(lon=-119, lat=34.5, period=1, radius=15)
        """
        
        args = {**self.args, **kwargs}
        sites = self.sites
        
        left, right = sites['lon'].min(), sites['lon'].max()
        bot, top = sites['lat'].min(), sites['lat'].max()
        ext = 0.3
        gmt_cmd = (f'range=-R{left-ext}/{right+ext}/{bot-ext}/{top+ext}\n'
                   f'model_region=-R{left}/{right}/{bot}/{top}\n'
                   f'gmt begin {str(self.directory)}/sites_location png,pdf A+s16c,H3,E600\n'
                   f'gmt set MAP_FRAME_TYPE plain\n'
                   f'gmt set FORMAT_FLOAT_OUT %.12lg\n'
                   f'gmt basemap -JM10c $range -Baf -BSWen\n'
                   f'gmt grdimage -C @earth_relief_{args["resolution"]}s -I+d\n'
                   f'gmt colorbar -DjMR+o-1c/0+m -I0.3 -G0/NaN -Bx+l"Elevation (m)"\n'
                   f'gmt plot sites_location.txt -St0.4c -Gcyan -W0.2p,black\n'
                   f'echo {self.lon} {self.lat} | gmt plot -Sa0.6c -Gred\n'
                   f'gmt coast -JM10c -Swhite\n'
                   f'# gmt plot borders.txt -W3p,lightblue -L\n'
                   f'# gmt plot -W2p,black -L << EOF\n'
                   f'# {left} {bot}\n'
                   f'# {right} {bot}\n'
                   f'# {right} {top}\n'
                   f'# {left} {top}\n'
                   f'# EOF\n'
                   f'gmt end')
        with open('plot_gof.gmt', 'w') as fid_gmt:
            fid_gmt.write(gmt_cmd)
        cmd = ['bash', 'plot_gof.gmt']
        call_sub(cmd, shell=False, timeout=args['timeout'])
    

    def get_run_id(self, site_name):
        cmd = (
            f'select R.Run_ID '
            f'from CyberShake_Runs R, CyberShake_Sites S '
            f'where R.Site_ID=S.CS_Site_ID '
            f'and R.Study_ID={self.STUDY_ID} '
            f'and S.CS_Short_Name={repr(site_name)} '
        )
        run_id = pd.read_sql(cmd, self.conn)
        return run_id.iloc[0].item()


    def run(self, **kwargs):
        args = {**self.args, **kwargs}
        try:
            len(args['period'])
            period = args['period'][0]
        except:
            period = args['period']
    
        count = len(list(self.directory.glob('*Disagg*txt')))
        if not count or args['force_update']:
            for site_name in self.sites['name']:
                run_id = self.get_run_id(site_name)
                cmd = [
                    f'./disagg_plot_wrapper.sh',
                    f'--run-id {run_id}',
                    f'--component RotD50',
                    f'--period {period}',
                    f'--probs {args["prob"]} --erf-file MeanUCERF.xml --type txt',
                    f'-o {self.directory}',
                ]
                try:
                    print(call_sub(cmd, timeout=args['timeout']))
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
                    print(f"disagg_plot_wrapper.sh failed for {site_name}: ", e)
                    continue
                print(f'Site {site_name} queried')

        dfs = []
        df_sources = []
        sites = self.sites
        for f in self.directory.glob('*Disagg*txt'):
            file_name = str(f)
            site_name = file_name.split('/')[1].split('_')[0]
            if sites[sites['name']==site_name]['dist'].item() < args['radius']:
                x, y = self.get_disagg(file_name)
                dfs.append(x)
                df_sources.append(y)
        nsite = len(dfs)
        dfs = pd.concat(dfs)
        df_sources = pd.concat(df_sources)\
           .groupby(['source_desc', 'source_id', 'erf_id', 'rup_id',
                     'rup_var_id', 'vel_model_id'])[['Contribution', 'ExceedRate']]\
           .sum().reset_index().sort_values('Contribution', ascending=False)
        df_sources.astype(object).set_index('source_desc', inplace=True)


        fig = plt.figure(figsize=plt.figaspect(0.33), dpi=300)
        fig.subplots_adjust(wspace=0.5)
        fig.suptitle(f"Disaggregate {nsite} sites within {args['radius']:.2g}km centered"
                     f" at ({self.lon:.2f}, {self.lat:.2f}), based on PSA-{period:.2g}s")
        ax = fig.add_subplot(1, 3, 1, projection='3d')
        res = dfs.groupby(level=0).min()
        self.plot_disagg(res, title="Min", ax=ax)
        
        ax = fig.add_subplot(1, 3, 2, projection='3d')
        res = dfs.groupby(level=0).mean()
        self.plot_disagg(res, title="Mean", ax=ax)
        
        ax = fig.add_subplot(1, 3, 3, projection='3d')
        res = dfs.groupby(level=0).max()
        ax.tick_params(axis='both', which='major', pad=0)
        self.plot_disagg(res, title="Max", ax=ax)

        #fig, ax = plt.subplots(1, len(colors), dpi=500, gridspec_kw={'wspace': 0.8})
        colors = mpl.cm.bwr_r(np.linspace(0, 1, len(res.columns) - 3))
        x, y = np.linspace(0, 1, 40), np.linspace(0, 1, 40)
        z = np.ones((40, 40))
        for i, color in enumerate(colors):
            ax = fig.add_subplot(position=[0.1 + 0.8 * i/len(colors), -0.06, 0.8/len(colors), 0.06])
            ax.pcolormesh(x, y, z, color=color, shading='auto')
            ax.set_axis_off()
            ax.set_aspect(1/2)
            pos = list(ax.get_position().bounds)
            x_text = pos[0] #+ pos[2] / 3 
            y_text = pos[1] - 1.2
            ax.text(0.5, y_text, res.columns[i + 2], ha='center', fontsize=12)
        fig.savefig(Path(self.directory, 'disagg.png'), dpi=600, bbox_inches='tight', pad_inches=0.05)
        df_sources.to_csv(Path(self.directory, 'scenarios.csv'), float_format='%.3g', index=False)
        print(f'\nThe {args["topk"]} most severe scenarios are:\n')
        print(f'{df_sources[:args["topk"]].set_index("source_desc")}')
        return dfs, df_sources
    
    # def plot_SA_map(self, **kwargs):
    #     args = {**self.args, **kwargs}
    #     if any('spatial' in k for k in args.keys()):
    #         cmd = (
    #             f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #             f"--study {args['study']} --period {args['period']} ",
    #             f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
    #             f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} ",
    #             f"--spatial-corr-fields {args['spatial_corr_fields']} ",
    #             f"--spacing {args['spacing']} -o {self.directory}",
    #             )
    #     else:
    #         cmd = (
    #             f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #             f"--study {args['study']} --period {args['period']} ",
    #             f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
    #             f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} ",
    #             f"--spacing {args['spacing']} -o {self.directory}",
    #             )
    #     print("".join(cmd))
    #     call_sub("".join(cmd), shell=True, timeout=args['timeout'])

    # def plot_SA_map(self, **kwargs):
    #     args = {**self.args, **kwargs}
    #     if any('spatial' in k for k in args.keys()):
    #         # cmd = (
    #         #     f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #         #     f"--study {args['study']} --period {args['period']} ",
    #         #     f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
    #         #     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} ",
    #         #     f"--spatial-corr-fields {args['spatial_corr_fields']} ",
    #         #     f"--spacing {args['spacing']} -o {self.directory}",
    #         #     )
    #         # Raw Cybershake
    #         cmd = (
    #             f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #             f"--study {args['study']} --period {args['period']} ",
    #             f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
    #             f"--rupture-var-id {args['rup_var_id']} ",
                
    #             f"--spacing {args['spacing']} -o {self.directory}",
    #             )
            
    #         # cmd = (
    #         #     f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #         #     f"--study {args['study']} --period {args['period']} ",
    #         #     f"--source-id {243} --rupture-id {8} ",
    #         #     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} ",
    #         #     f"--spacing {args['spacing']} -o {self.directory}",
    #         #     )


    #     else:
    #         cmd = (
    #             f"java -Xmx2G -cp opensha-cybershake-all.jar ",
    #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
    #             f"--study {args['study']} --period {args['period']} ",
    #             f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
    #             f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} ",
    #             f"--spacing {args['spacing']} -o {self.directory}",
    #             )
            
    #     print("".join(cmd))
    #     call_sub("".join(cmd), shell=True, timeout=args['timeout'])

         

       



# #Working section                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
#     def plot_SA_map(self, **kwargs):
#         args = {**self.args, **kwargs}
#         if any('spatial' in k for k in args.keys()):
#             rup_ID_try = 8
#             while rup_ID_try == 8:
#                 try:
#                     # cmd = (
#                     #     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     #     f"--study {args['study']} --period {args['period']} "
#                     #     f"--source-id {args['source_id']}  --rupture-id {args['rup_id']} "
#                     #     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                     #     f"--spacing {args['spacing']} -o {self.directory}"
#                     #  )
#                 #gmpe at cs sites
#                     # cmd = (
#                     #     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     #     f"--study {args['study']} --period {args['period']} "
#                     #     f"--source-id {args['source_id']}  --rupture-id {args['rup_id']} "
#                     #     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                     #     f"--calc-gmpe-at-cs-sites "
#                     #     f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                     #     f"--spacing {args['spacing']} -o {self.directory}"
#                     # )
#                 # Raw cybershake
#                     cmd = (
#                         f"java -Xmx2G -cp opensha-cybershake-all.jar ",
#                         f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
#                         f"--study {args['study']} --period {args['period']} ",
#                         f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
#                         f"--rupture-var-id {args['rup_var_id']} ",
#                         f"--spatial-corr-fields {args['spatial_corr_fields']} ",
#                         # f"--colorbar-min {0} ",
#                         # f"--colorbar-max {0.6} ",
#                         f"--spatial-corr-debug ",  
#                         f"--download-interpolated ", 
#                         f"--spacing {args['spacing']} -o {self.directory}",
#                         )
#                 # Interp grid
#                     # cmd = (
#                     #     f"java -Xmx2G -cp opensha-cybershake-all.jar ",
#                     #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
#                     #     f"--study {args['study']} --period {args['period']} ",
#                     #     f"--source-id {args['source_id']} --rupture-id {args['rup_id']} ",
#                     #     f"--download-interpolated ",
#                     #     f"--rupture-var-id {args['rup_var_id']} ",
#                     #     f"--spacing {args['spacing']} -o {self.directory} ", 
#                     #     )
                   
                    
#                     print("".join(cmd))
#                     call_sub("".join(cmd), shell=True, timeout=args['timeout'])
#                     break
#                 except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                     print(f"CyberShakeScenarioShakeMapGenerator failed: {e}")

#                     try:
#                         cmd = (
#                             f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                             f"--study {args['study']} --period {args['period']} "
#                             f"--source-id {args['source_id']} --rupture-id {4} "
#                             f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                             f"--spacing {args['spacing']} -o {self.directory}"
#                         )
#                     # gmpe at cs sites
#                         # cmd = (
#                         #     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                         #     f"--study {args['study']} --period {args['period']} "
#                         #     f"--source-id {args['source_id']}  --rupture-id {4} "
#                         #     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                         #     f"--calc-gmpe-at-cs-sites "
#                         #     f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                         #     f"--spacing {args['spacing']} -o {self.directory}"
#                         # )
#                     # raw cybershake
#                         # cmd = (
#                         # f"java -Xmx2G -cp opensha-cybershake-all.jar ",
#                         # f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
#                         # f"--study {args['study']} --period {args['period']} ",
#                         # f"--source-id {args['source_id']} --rupture-id {4} ",
#                         # f"--rupture-var-id {args['rup_var_id']} ",
#                         # f"--spatial-corr-fields {args['spatial_corr_fields']} ", 
#                         # f"--spacing {args['spacing']} -o {self.directory}",
#                         # )
#                     # Interp grid
#                         # cmd = (
#                         #     f"java -Xmx2G -cp opensha-cybershake-all.jar ",
#                         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator ",
#                         #     f"--study {args['study']} --period {args['period']} ",
#                         #     f"--source-id {args['source_id']} --rupture-id {4} ",
#                         #     f"--download-interpolated ",
#                         #     f"--rupture-var-id {args['rup_var_id']} ",
#                         #     f"--spacing {args['spacing']} -o {self.directory} ", 
#                         #     )

                        
#                         print("".join(cmd))
#                         call_sub("".join(cmd), shell=True, timeout=args['timeout'])
#                         break
#                     except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                         print(f"CyberShakeScenarioShakeMapGenerator retry failed: {e}")
#                         continue
#         else:
#             try:
#                 cmd = (
#                     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     f"--study {args['study']} --period {args['period']} "
#                     f"--source-id {args['source_id']} --rupture-id {args['rup_id']} "
#                     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                     f"--spacing {args['spacing']} -o {self.directory}"
#                 )

                
#                 print("".join(cmd))
#                 call_sub("".join(cmd), shell=True, timeout=args['timeout'])
#             except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                 print(f"CyberShakeScenarioShakeMapGenerator failed: {e}")
    

#    #Rupture Varations 
#     def plot_SA_map(self, **kwargs):
#         args = {**self.args, **kwargs}
#         if any('spatial' in k for k in args.keys()):
#             rup_ID_try = 8
#             while rup_ID_try == 8:
#                 try:
                    
#                     # First attempt with rupture ID 8
#                     for i in range(20,25):
#                         # cmd = (
#                         #     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                         #     f"--study {args['study']} --period {args['period']} "
#                         #     f"--source-id {args['source_id']} --rupture-id {rup_ID_try} "
#                         #     f"--download-interpolated "
#                         #     f"--rupture-var-id {i} "
#                         #     f"--spatial-corr-debug ",
#                         #     f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                         #     f"--spacing {args['spacing']} -o {self.directory}"
#                         # )
#                         cmd = (
#                             f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                             f"--study {args['study']} --period {args['period']} "
#                             f"--source-id {223} --rupture-id {3} "
#                             f"--download-interpolated "
#                             f"--rupture-var-id {i} "
#                             f"--spatial-corr-debug ",
#                             f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                             f"--spacing {args['spacing']} -o {self.directory}"
#                         )
                        
#                         print("".join(cmd))
#                         call_sub("".join(cmd), shell=True, timeout=args['timeout'])
                    
#                     # Break loop if successful
#                     break
#                 except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                     print(f"CyberShakeScenarioShakeMapGenerator failed with rupture ID {rup_ID_try}: {e}")
#                     try:
#                         sys.exit(0)
#                     # try:
#                     #     # Retry with rupture ID 4 if initial attempt fails
#                     #     for i in range(45):
#                     #         cmd = (
#                     #             f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     #             f"--study {args['study']} --period {args['period']} "
#                     #             f"--source-id {args['source_id']} --rupture-id {4} "
#                     #             f"--download-interpolated "
#                     #             f"--rupture-var-id {i} "
#                     #             f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                     #             f"--spacing {args['spacing']} -o {self.directory}"
#                     #         )
                            
#                     #         print("".join(cmd))
#                     #         call_sub("".join(cmd), shell=True, timeout=args['timeout'])
                        
#                         # Break loop if successful
#                         break
#                     except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                         print(f"CyberShakeScenarioShakeMapGenerator retry failed with rupture ID 4: {e}")
#                         continue
#         else:
#             try:
#                 cmd = (
#                     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     f"--study {args['study']} --period {args['period']} "
#                     f"--source-id {args['source_id']} --rupture-id {args['rup_id']} "
#                     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                     f"--spacing {args['spacing']} -o {self.directory}"
#                 )
#                 print("".join(cmd))
#                 call_sub("".join(cmd), shell=True, timeout=args['timeout'])
#             except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                 print(f"CyberShakeScenarioShakeMapGenerator failed: {e}")

#    #Northridge 
#     def plot_SA_map(self, **kwargs):
#         args = {**self.args, **kwargs}
#         if any('spatial' in k for k in args.keys()):
#             rup_ID_try = 8
#             while rup_ID_try == 8:
#                 try:
                    
#                     # First attempt with rupture ID 8
#                     for i in range(63):
#                         # cmd = (
#                         #     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                         #     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                         #     f"--study {args['study']} --period {args['period']} "
#                         #     f"--source-id {args['source_id']} --rupture-id {rup_ID_try} "
#                         #     f"--download-interpolated "
#                         #     f"--rupture-var-id {i} "
#                         #     f"--spatial-corr-debug ",
#                         #     f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                         #     f"--spacing {args['spacing']} -o {self.directory}"
#                         # )
#                         cmd = (
#                             f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                             f"--study {args['study']} --period {args['period']} "
#                             f"--source-id {0} --rupture-id {0} "
#                             f"--rupture-var-id {0} "
#                             f"--download-interpolated "
#                             f"--spatial-corr-debug "
#                             f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                             f"--spacing {args['spacing']} -o {self.directory}"
#                         )
                        
#                         print("".join(cmd))
#                         call_sub("".join(cmd), shell=True, timeout=args['timeout'])
                    
#                     # Break loop if successful
#                     break
#                 except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                     print(f"CyberShakeScenarioShakeMapGenerator failed with rupture ID {rup_ID_try}: {e}")
#                     try:
#                         sys.exit(0)
#                     # try:
#                     #     # Retry with rupture ID 4 if initial attempt fails
#                     #     for i in range(45):
#                     #         cmd = (
#                     #             f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     #             f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     #             f"--study {args['study']} --period {args['period']} "
#                     #             f"--source-id {args['source_id']} --rupture-id {4} "
#                     #             f"--download-interpolated "
#                     #             f"--rupture-var-id {i} "
#                     #             f"--spatial-corr-fields {args['spatial_corr_fields']} "
#                     #             f"--spacing {args['spacing']} -o {self.directory}"
#                     #         )
                            
#                     #         print("".join(cmd))
#                     #         call_sub("".join(cmd), shell=True, timeout=args['timeout'])
                        
#                         # Break loop if successful
#                         break
#                     except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                         print(f"CyberShakeScenarioShakeMapGenerator retry failed with rupture ID 4: {e}")
#                         continue
#         else:
#             try:
#                 cmd = (
#                     f"java -Xmx2G -cp opensha-cybershake-all.jar "
#                     f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
#                     f"--study {args['study']} --period {args['period']} "
#                     f"--source-id {args['source_id']} --rupture-id {args['rup_id']} "
#                     f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
#                     f"--spacing {args['spacing']} -o {self.directory}"
#                 )
#                 print("".join(cmd))
#                 call_sub("".join(cmd), shell=True, timeout=args['timeout'])
#             except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
#                 print(f"CyberShakeScenarioShakeMapGenerator failed: {e}")


    def plot_SA_map(self, **kwargs):
        args = {**self.args, **kwargs}
        if any('spatial' in k for k in args.keys()):
            for i in range(45):
                try:
                    # Command for generating maps
                    cmd = (
                        f"java -Xmx2G -cp opensha-cybershake-all.jar "
                        f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
                        f"--study {args['study']} --period {args['period']} "
                        f"--source-id {args['source_id']} --rupture-id {args['rup_id']} "
                        f"--download-interpolated "
                        f"--rupture-var-id {i} "
                        f"--spatial-corr-fields {args['spatial_corr_fields']} "
                        f"--spacing {args['spacing']} -o {self.directory}"
                    )
                    
                    print("Executing command:", cmd)
                    subprocess.run(cmd, shell=True, check=True, timeout=args.get('timeout', 60))

                except subprocess.CalledProcessError as e:
                    print(f"CyberShakeScenarioShakeMapGenerator failed for rupture-var-id {i}: {e}")
                    # Handle specific errors if necessary or skip to the next iteration

                except subprocess.TimeoutExpired as e:
                    print(f"CyberShakeScenarioShakeMapGenerator timed out for rupture-var-id {i}: {e}")
                    # Handle timeout specifically if necessary or skip to the next iteration

        else:
            try:
                cmd = (
                    f"java -Xmx2G -cp opensha-cybershake-all.jar "
                    f"org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator "
                    f"--study {args['study']} --period {args['period']} "
                    f"--source-id {args['source_id']} --rupture-id {args['rup_id']} "
                    f"--rupture-var-id {args['rup_var_id']} --gmpe {args['gmpe']} "
                    f"--spacing {args['spacing']} -o {self.directory}"
                )
                print("Executing command:", cmd)
                subprocess.run(cmd, shell=True, check=True, timeout=args.get('timeout', 60))

            except subprocess.CalledProcessError as e:
                print(f"CyberShakeScenarioShakeMapGenerator failed: {e}")
            
            except subprocess.TimeoutExpired as e:
                print(f"CyberShakeScenarioShakeMapGenerator timed out: {e}")



if __name__ == '__main__':
    if not Path('disagg_plot_wrapper.sh').exists() \
            or not Path('run_opensha.sh').exists(): 
        with open('disagg_plot_wrapper.sh', 'w') as fid:
            cmd = """
#!/bin/bash

set -o errexit

if [[ $# -lt 1 ]];then
	echo "must supply class to run (and optional arguments)"
	exit 2
fi

java=`which java`

classpath="opensha-cybershake-all.jar"
#$java -Xms512M -Xmx2G -cp $classpath -Dcybershake.db.host=updated_study_15_12.sqlite -Xmx12g org.opensha.sha.cybershake.plot.DisaggregationPlotter $@
$java -Xms512M -Xmx2G -cp $classpath -Xmx12g org.opensha.sha.cybershake.plot.DisaggregationPlotter $@
"""
    if not Path('opensha-cybershake-all.jar').exists():
        print("Jar file (opensha-cybershake-all.jar) not found!\nAborting!")
        import sys 
        sys.exit(-1)
    disagg = Disagg()

    

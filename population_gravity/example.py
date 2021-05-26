import numpy as np
from population_gravity import Model


class TestProjectedOutputs:
    """Test configuration integrity."""

    def __init__(self, au, ar, bu, br, kd, rid, state_name, ssp, root_dir, hist_year, start_year, end_year,
                 urban_pop_n=None, rural_pop_n=None):

        self.au = au
        self.ar = ar
        self.bu = bu
        self.br = br
        self.kd = kd
        self.rid = rid
        self.hist_year = hist_year
        self.start_year = start_year
        self.end_year = end_year
        self.urban_pop_n = urban_pop_n
        self.rural_pop_n = rural_pop_n

        if hist_year <= 2010:
            self.HIST_RURAL_RASTER = f'{root_dir}/inputs/{state_name}_rural_{self.hist_year}_1km.tif'
            self.HIST_URBAN_RASTER = f'{root_dir}/inputs/{state_name}_urban_{self.hist_year}_1km.tif'

        else:
            self.HIST_RURAL_RASTER = f'{root_dir}/inputs/{state_name}_{ssp}_rural_{self.hist_year}_1km.tif'
            self.HIST_URBAN_RASTER = f'{root_dir}/inputs/{state_name}_{ssp}_urban_{self.hist_year}_1km.tif'

        self.HIST_SUITABILITY = f'{root_dir}/inputs/{state_name}_mask_short_term.tif'
        self.GRID_COORD_FILE = f'{root_dir}/inputs/{state_name}_coordinates.csv'
        self.PROJ_POP_FILE = f'{root_dir}/inputs/{state_name}_{ssp}_popproj.csv'

        self.ONE_D_IND_FILE = f'{root_dir}/inputs/{state_name}_within_indices.txt'
        self.OUTPUT_DIRECTORY = f'{root_dir}/outputs'

        self.SCENARIO = ssp
        self.STATE_NAME = state_name

    def run_array1d_outputs(self):
        """Test for 1D array outputs."""

        run = Model(grid_coordinates_file=self.GRID_COORD_FILE,
                    base_rural_pop_raster=self.HIST_RURAL_RASTER,
                    base_urban_pop_raster=self.HIST_URBAN_RASTER,
                    historical_suitability_raster=self.HIST_SUITABILITY,
                    projected_population_file=self.PROJ_POP_FILE,
                    urban_pop_proj_n=self.urban_pop_n,
                    rural_pop_proj_n=self.rural_pop_n,
                    one_dimension_indices_file=self.ONE_D_IND_FILE,
                    output_directory=self.OUTPUT_DIRECTORY,
                    alpha_urban=self.au,
                    alpha_rural=self.ar,
                    beta_urban=self.bu,
                    beta_rural=self.br,
                    kernel_distance_meters=self.kd,
                    scenario=self.SCENARIO,
                    state_name=self.STATE_NAME,
                    historic_base_year=self.hist_year,
                    projection_year=self.start_year,
                    projection_end_year=self.end_year,
                    time_step=10,
                    write_raster=True,
                    write_array1d=True,
                    write_array2d=False,
                    write_csv=False,
                    write_suitability=False,
                    write_logfile=False,
                    run_number=self.rid,
                    output_total=False)

        run.downscale()


if __name__ == '__main__':

    import os
    import pandas as pd

    ssp = 'ssp2'
    root_dir = '/Users/d3y010/projects/population/california'
    # root_dir = '/Users/d3y010/projects/population/zoraghein-oneill_population_gravity_inputs_outputs'
    # root_dir = '/Users/d3y010/Desktop/staging'
    input_dir = os.path.join(root_dir, 'inputs')
    # input_dir = root_dir

    historic_year = 2000
    start_year = 2010
    end_year = 2010

    # use existing params from publication?
    use_params = False

    # number of samples if running batch
    n_samples = 1000

    # run id
    run_id = ''

    # run array
    run_array = range(0, n_samples)

    # get a list of all states
    # neighbors_file = '/Users/d3y010/repos/github/population_gravity/population_gravity/data/neighboring_states_150km.csv'
    # ndf = pd.read_csv(neighbors_file)
    # state_list = ndf['target_state'].unique().tolist()

    # already ran
    # ran_states = [i.split('_1km')[0] for i in os.listdir(os.path.join(root_dir, 'outputs')) if '_1km_ssp2_urban_' in i]
    # state_list = [i for i in state_list if i not in ran_states]

    state_list = ['california']

    # MANUAL OVERIDE OF POPULATION PROJECTION VALUES; SET TO None if using file
    # urban_pop_n = None
    # rural_pop_n = None
    urban_pop_n = 35221640
    rural_pop_n = 1872743

    for state_name in state_list:
        print(f"Processing:  {state_name}")

        # root_dir = f'/Users/d3y010/projects/population/zoraghein-oneill_population_gravity_inputs_outputs/{state_name}'
        # input_dir = os.path.join(root_dir, 'inputs')

        if use_params:

            # pdf = pd.read_csv(os.path.join(input_dir, f"{state_name}_{ssp}_params.csv"))
            pdf = pd.read_csv(f'/Users/d3y010/Desktop/validation_params/{state_name}_ssp2_params_{historic_year}to{start_year}.csv')

            alpha_urban = pdf['Alpha_Urban'].values[0]
            alpha_rural = pdf['Alpha_Rural'].values[0]
            beta_urban = pdf['Beta_Urban'].values[0]
            beta_rural = pdf['Beta_Rural'].values[0]
            kernel_distance_meters = 100000

            pop = TestProjectedOutputs(au=alpha_urban,
                                       ar=alpha_rural,
                                       bu=beta_urban,
                                       br=beta_rural,
                                       kd=kernel_distance_meters,
                                       rid=run_id,
                                       state_name=state_name,
                                       ssp=ssp,
                                       root_dir=root_dir,
                                       hist_year=historic_year,
                                       start_year=start_year,
                                       end_year=end_year,
                                       urban_pop_n=urban_pop_n,
                                       rural_pop_n=rural_pop_n)

            pop.run_array1d_outputs()

        else:

            params_arr = os.path.join(root_dir, 'outputs', 'lhs', f"lhs_{n_samples}_sample.npy")
            arr_params = np.load(params_arr)

            df_params = pd.DataFrame({'alpha_urban': arr_params[:, 0],
                                      'alpha_rural': arr_params[:, 1],
                                      'beta_urban': arr_params[:, 2],
                                      'beta_rural': arr_params[:, 3],
                                      'kernel_distance_meters': arr_params[:, 4]})

            df_params['n'] = df_params.index

            for run_id in range(n_samples):

                x = df_params.loc[df_params['n'] == run_id]

                alpha_urban = x['alpha_urban'].values[0]
                alpha_rural = x['alpha_rural'].values[0]
                beta_urban = x['beta_urban'].values[0]
                beta_rural = x['beta_rural'].values[0]
                kernel_distance_meters = x['kernel_distance_meters'].values[0]

                pop = TestProjectedOutputs(au=alpha_urban,
                                           ar=alpha_rural,
                                           bu=beta_urban,
                                           br=beta_rural,
                                           kd=kernel_distance_meters,
                                           rid=run_id,
                                           state_name=state_name,
                                           ssp=ssp,
                                           root_dir=root_dir,
                                           hist_year=historic_year,
                                           start_year=start_year,
                                           end_year=end_year,
                                           urban_pop_n=urban_pop_n,
                                           rural_pop_n=rural_pop_n)

                pop.run_array1d_outputs()

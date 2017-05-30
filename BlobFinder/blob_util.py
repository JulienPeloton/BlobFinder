import numpy as np

class normalise_parser(object):
    """
    Generic class to handle parsered file.
    Mostly contain basic static functions.
    """
    def __init__(self, config_dict):
        """
        Parameters
        ----------
            * config_dict: dic, dictionary coming from the ini file.
        """
        self.config_dict = config_dict

    @staticmethod
    def floatise_it(entry):
        return float(entry)

    @staticmethod
    def intise_it(entry):
        return int(entry)

    @staticmethod
    def boolise_it(dic, entry):
        if entry in dic:
            out = 'True' in dic[entry]
        else:
            out = False
        return out

    @staticmethod
    def normalise_array(array, func):
        array_out = np.array([func(i) for i in array.split()])
        if len(array_out) > 0:
            return array_out
        else:
            return None

    @staticmethod
    def make_dic(keys, array, func):
        values = np.array([func(i) for i in array.split()])
        return {i: j for i, j in zip(keys, values)}

class normalise_instrument_parser(normalise_parser):
    """
    Class to handle foreground parser.
    It converts initial dictionary into an object.
    """
    def __init__(self, config_dict):
        """
        Parameters
        ----------
            * config_dict: dic, dictionary coming from the ini file.
        """
        normalise_parser.__init__(self, config_dict)

        ## Names
        self.output_path = config_dict['output_path']
        self.fiducial_cmb = config_dict['fiducial_cmb']
        ## Booleans

        ## Integers
        self.npix = self.intise_it(config_dict['npix'])
        self.Number_of_Sources = self.intise_it(
            config_dict['number_of_sources'])
        self.Number_of_Sources_EX = self.intise_it(
            config_dict['number_of_sources_ex'])
        self.Number_of_SZ_Clusters = self.intise_it(
            config_dict['number_of_sz_clusters'])
        self.n_iterations = self.intise_it(config_dict['n_iterations'])
        self.seed = self.intise_it(config_dict['seed'])

        ## float
        self.pix_size = self.floatise_it(config_dict['pix_size'])
        self.snr_threshold = self.floatise_it(config_dict['snr_threshold'])
        self.white_noise_level = self.floatise_it(
            config_dict['white_noise_level'])
        self.atmospheric_noise_level = self.floatise_it(
            config_dict['atmospheric_noise_level'])
        self.one_over_f_noise_level = self.floatise_it(
            config_dict['one_over_f_noise_level'])
        self.beam_size_fwhm = self.floatise_it(config_dict['beam_size_fwhm'])
        self.Amplitude_of_Sources = self.floatise_it(
            config_dict['amplitude_of_sources'])
        self.Amplitude_of_Sources_EX = self.floatise_it(
            config_dict['amplitude_of_sources_ex'])
        self.Mean_Amplitude_of_SZ_Clusters = self.floatise_it(
            config_dict['mean_amplitude_of_sz_clusters'])
        self.SZ_beta = self.floatise_it(config_dict['sz_beta'])
        self.SZ_Theta_core = self.floatise_it(config_dict['sz_theta_core'])
        self.X_width = self.npix * self.pix_size / 60.
        self.Y_width = self.npix * self.pix_size / 60.

import argparse
import ConfigParser
import os
import sys

import numpy as np
import pylab as pl

import blob_util
import blob

def addargs(parser):
    ''' Parse command line arguments '''
    parser.add_argument(
        '-setup_instrument', dest='setup_instrument',
        required=True,
        help='Configuration file with parameters.')
    parser.add_argument(
        '--plot', dest='plot',
        action='store_true',
        help='Plot output')

def grabargs(args_param=None):
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser(
        description='Package to find point sources in sky maps.')
    addargs(parser)
    args = parser.parse_args(args_param)

    Config = ConfigParser.ConfigParser()
    Config.read(args.setup_instrument)
    environment = blob_util.normalise_instrument_parser(
        Config._sections['SimulationParameters'])

    ## Create output folder
    if not os.path.exists(environment.output_path):
        os.makedirs(environment.output_path)

    return args, environment

if __name__ == '__main__':
    args_param = None
    args, p = grabargs(args_param)

    # ###############################################
    # ## Create a CMB map
    # ###############################################
    ell, DlTT = np.loadtxt(p.fiducial_cmb, usecols=(0, 1), unpack=True)
    CMB_T = blob.make_CMB_T_map(p.npix, p.pix_size, ell, DlTT, seed=p.seed)

    total_map_plus_noise_array = []
    SN_map_array = []
    pixels_array = []
    for CMB in [CMB_T]:
        PSMap = blob.Poisson_source_component(
            p.npix, p.pix_size, p.Number_of_Sources,
            p.Amplitude_of_Sources, seed=p.seed)

        PSMap += blob.Exponential_source_component(
            p.npix, p.pix_size, p.Number_of_Sources_EX,
            p.Amplitude_of_Sources_EX, seed=p.seed)

        ## make an SZ map
        SZMap, SZCat = blob.SZ_source_component(
            p.npix, p.pix_size, p.Number_of_SZ_Clusters,
            p.Mean_Amplitude_of_SZ_Clusters, p.SZ_beta,
            p.SZ_Theta_core, False, seed=p.seed)

        ###############################################
        ## Add all components to get the sky map at a single frequency
        ###############################################
        total_map = CMB + PSMap + SZMap

        ###############################################
        ## Convolved with the beam of the instrument
        ###############################################
        CMB_T_convolved = blob.convolve_map_with_gaussian_beam(
            p.npix, p.pix_size, p.beam_size_fwhm, total_map)

        ###############################################
        ## Add noise
        ###############################################
        Noise = blob.make_noise_map(
            p.npix, p.pix_size, p.white_noise_level,
            p.atmospheric_noise_level, p.one_over_f_noise_level, seed=p.seed)

        total_map_plus_noise = CMB_T_convolved + Noise
        total_map_plus_noise_array.append(total_map_plus_noise)

        FT_noise_covar = blob.construct_noise_covar(
            p.npix, p.pix_size, p.beam_size_fwhm,
            ell, DlTT,
            p.Number_of_Sources, p.Amplitude_of_Sources,
            p.Number_of_Sources_EX, p.Amplitude_of_Sources_EX,
            p.Number_of_SZ_Clusters, p.Mean_Amplitude_of_SZ_Clusters,
            p.SZ_beta, p.SZ_Theta_core,
            p.white_noise_level, p.atmospheric_noise_level,
            p.one_over_f_noise_level, p.n_iterations, seed=57467536)

        ## construct the beam and cluster profile for
        ## the numerator of the matched filter
        ## this is the filtering we did on the map
        beam_and_filt = blob.make_2d_gaussian_beam(
            p.npix, p.pix_size, p.beam_size_fwhm)

        ## this is the signal we are looking for
        cluster_profile = blob.beta_function(
            p.npix, p.pix_size, p.SZ_beta, p.SZ_Theta_core)

        ## Apply the matched filter to our map
        filtered_map = blob.matched_filter(
            total_map_plus_noise, beam_and_filt,
            cluster_profile, FT_noise_covar)

        ## make a S/N map
        SN_map = filtered_map / np.std(filtered_map)
        SN_map_array.append(SN_map)

        pixels = np.where(SN_map > p.snr_threshold)
        pixels_array.append(pixels)
        print ' '
        print "%d sources with SNR greater than %.2f: " % (
            len(pixels[0]), p.snr_threshold), pixels
        true_sources = blob.sort_sources(pixels, threshold=5)

    if args.plot:
        ## make a few plots
        fig, ax = pl.subplots(1, 2, figsize=(13, 5))
        blob.plot_sky_map(
            total_map_plus_noise_array[0],
            -400, 400,
            p.X_width, p.Y_width,
            ax[0], title='I Sky+beam+noise')

        blob.plot_snr_map(
            SN_map_array[0],
            -5, 10,
            p.X_width, p.Y_width,
            ax=ax[1], title='Sky+beam+noise (SNR)',
            sources_position=true_sources)
        pl.show()

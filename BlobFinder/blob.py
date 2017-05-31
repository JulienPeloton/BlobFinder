from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import pylab as pl
import numpy as np


def sort_sources(pixels, threshold=5):
    """
    Remove pixels which belong to the same source.

    Parameters
    -----------
        * pixels: 2D array, x and y coordinates of the pixels having high SNR
        * threshold: int, the minimum separation length (in pixel) to
            consider two pixels as belonging to different sources.
    """
    true_sources_x = []
    true_sources_y = []
    for pos, (x, y) in enumerate(zip(pixels[0], pixels[1])):
        if pos == 0:
            true_sources_x.append(x)
            true_sources_y.append(y)
            continue
        if np.any(np.abs(true_sources_x - x) < threshold) and \
        np.any(np.abs(true_sources_y - y) < threshold):
            continue
        else:
            true_sources_x.append(x)
            true_sources_y.append(y)
    return np.array([true_sources_x, true_sources_y])

def matched_filter(input_map, beam_and_filt, signal_profile, FT_noise_covar):
    """
    Parameters
    -----------
        * input_map: 2D array, the map we are processing
        * beam_and_filt: 2D array, the beam convolved with
            any map filtering, in real space.
        * signal_profile: 2D array, the shape of
            the signal we are looking for, in real spcae
        * FT_noise_covar: 2D array, the B_N_{ap}^2 + N_{ins}^2 in fourier space
    """
    ## tranform beam_and_filt to fourier space
    FT_beam_and_filt = np.fft.fft2(np.fft.fftshift(beam_and_filt))

    ## tranform cluster_profile to fourier space
    FT_signal = np.fft.fft2(np.fft.fftshift(signal_profile))

    ## define the matched filter function
    psi = FT_beam_and_filt * FT_signal / FT_noise_covar

    ## filter the map
    filtered = psi * np.fft.fft2(np.fft.fftshift(input_map))
    ## center the filter
    filtered = np.fft.fftshift(np.fft.ifft2(filtered))
    ## change the data type to real
    filtered = np.real(filtered)

    return filtered

def plot_snr_map(
    Map_to_Plot, c_min, c_max, X_width, Y_width,
    ax=None, title='', sources_position=None):
    """
    Plot SNR map

    Parameters
    -----------
        * Map_to_Plot: 2D array, the map to plot
        * c_min/c_max: int, min/max for colorbar [SNR]
        * X/Y_width: int, horizontal/vertical size [pixel]
        * sources_position: 2D array, x and y coordinates of
            the pixels having high SNR
    """
    Npix = len(Map_to_Plot)

    im = ax.imshow(Map_to_Plot, interpolation='bilinear',
        origin='lower',cmap=pl.cm.viridis)

    cols = []
    rows = []
    if sources_position is not None:
        print '+---------------------------+'
        for col,row in zip(sources_position[0], sources_position[1]):
            if col in cols or row in rows:
                continue
            cols.append(col)
            rows.append(row)
            ax.scatter(row * X_width / Npix, col * Y_width / Npix, s=500,
                facecolors='none', edgecolors='black')
            print "SNR at [%.2f, %.2f] deg: " % (
                row * X_width / Npix,
                col * Y_width / Npix), Map_to_Plot[col][row]

    ax.set_title(title)
    im.set_clim(c_min, c_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = pl.colorbar(im, cax=cax)
    im.set_extent([0, X_width, 0, Y_width])
    ax.set_ylabel('$\Delta$Dec $[^\circ]$')
    ax.set_xlabel('$\Delta$RA $[^\circ]$')
    # cbar.set_label('matched filtered data [S/N]', rotation=270)

def convolve_map_with_gaussian_beam(N, pix_size, beam_size_fwhm, Map):
    """
    convolves a map with a Gaussian beam pattern.

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size [arcmin]
        * beam_size_fwhm: float, FWHM of the instrument [arcmin]
        * Map: 2D array, the map to be convolved
    """

    # make a 2d gaussian
    gaussian = make_2d_gaussian_beam(N, pix_size, beam_size_fwhm)

    # do the convolution
    FT_gaussian = np.fft.fft2(np.fft.fftshift(gaussian))
    FT_Map = np.fft.fft2(np.fft.fftshift(Map))
    convolved_map = np.fft.fftshift(np.real(np.fft.ifft2(FT_gaussian * FT_Map)))

    return convolved_map
  ###############################

def make_2d_gaussian_beam(N, pix_size, beam_size_fwhm):
    """
    Construct a 2D Gaussian beam

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * beam_size_fwhm: float, FWHM of the instrument
    """
     # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N) + .5 - N / 2.) * pix_size
    X = np.outer(ones,inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)

    # make a 2d gaussian
    beam_sigma = beam_size_fwhm / np.sqrt(8. * np.log(2))
    gaussian = np.exp(-.5 * (R / beam_sigma)**2.)
    gaussian = gaussian / np.sum(gaussian)

    return gaussian
  ###############################

def make_CMB_T_map(N,pix_size,ell,DlTT,seed=59843757):
    """
    Makes a realization of a simulated CMB sky map (flat sky approximation)

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * ell: 1D array, multipoles (from 2 to lmax)
        * DlTT: 1D array, fiducial CMB temperature
            power spectrum (normalised by l*(l+1)/2pi)
    """
    state_initial = np.random.RandomState(seed)

    # convert Dl to Cl
    ClTT = DlTT * 2 * np.pi / (ell * (ell + 1.))
    ClTT[0] = 0.
    ClTT[1] = 0.

    # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N) + .5 - N / 2.) / (N - 1.)
    X = np.outer(ones, inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)

    # now make a 2d CMB power spectrum
    ell_scale_factor = 2. * np.pi / (pix_size / 60. * np.pi / 180.)
    ell2d = R * ell_scale_factor
    ClTT_expanded = np.zeros(int(ell2d.max() + 1))
    ClTT_expanded[0:(ClTT.size)] = ClTT
    CLTT2d = ClTT_expanded[ell2d.astype(int)]

    # now make a realization of the CMB
    # with the given power spectrum in fourier space
    ramdomn_array_for_T = np.fft.fft2(state_initial.normal(0, 1, (N, N)))
    FT_2d = np.sqrt(CLTT2d) * ramdomn_array_for_T
    ## make a plot of the 2d cmb simulated map in
    # fourier space, note the x and y axis labels need to be fixed
    CMB_T = np.fft.ifft2(
        np.fft.fftshift(FT_2d)) /(pix_size / 60. * np.pi / 180.)
    CMB_T = np.real(CMB_T)

    return CMB_T

def plot_sky_map(Map_to_Plot, c_min, c_max, X_width, Y_width, ax=None, title=''):
    """
    Plot sky map

    Parameters
    -----------
        * Map_to_Plot: 2D array, the map to plot
        * c_min/c_max: int, min/max for colorbar [uK]
        * X/Y_width: int, horizontal/vertical size [pixel]
        * sources_position: 2D array, x and y coordinates of
            the pixels having high SNR
    """
    im = ax.imshow(
        Map_to_Plot, interpolation='bilinear',
        origin='lower', cmap=pl.cm.viridis)

    ax.set_title(title)
    im.set_clim(c_min,c_max)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = pl.colorbar(im, cax=cax)
    im.set_extent([0, X_width, 0, Y_width])
    ax.set_ylabel('$\Delta$Dec $[^\circ]$')
    ax.set_xlabel('$\Delta$RA $[^\circ]$')
    cbar.set_label('Temperature [uK]', rotation=270)

def Poisson_source_component(
    N, pix_size, Number_of_Sources, Amplitude_of_Sources, seed=5439058):
    """
    Generate a faint distribution of sources with a poisson
    distribution of brightness.

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * Number_of_Sources: int, number of faint sources to add
        * Amplitude_of_Sources: int, parameter for the Poisson distribution

    """
    PSMap = np.zeros([N,N])
    state_initial = np.random.RandomState(seed)
    i = 0
    while (i < Number_of_Sources):
        pix_x = int(N*state_initial.rand())
        pix_y = int(N*state_initial.rand())
        PSMap[pix_x,pix_y] += state_initial.poisson(int(Amplitude_of_Sources))
        i = i + 1

    return(PSMap)
  ###############################

def Exponential_source_component(
    N, pix_size, Number_of_Sources_EX, Amplitude_of_Sources_EX, seed=487483):
    """
    Generate small number of very bright sources with and exponentially
    falling source count.

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * Number_of_Sources_EX: int, number of bright sources to add
        * Amplitude_of_Sources_EX: int, parameter for the exponential distribution

    """
    PSMap = np.zeros([N,N])
    state_initial = np.random.RandomState(seed)
    i = 0
    while (i < Number_of_Sources_EX):
        pix_x = int(N*state_initial.rand())
        pix_y = int(N*state_initial.rand())
        PSMap[pix_x,pix_y] += state_initial.exponential(Amplitude_of_Sources_EX)
        i = i + 1

    return(PSMap)
  ###############################

def SZ_source_component(
    N, pix_size, Number_of_SZ_Clusters, Mean_Amplitude_of_SZ_Clusters,
    SZ_beta, SZ_Theta_core, do_plots, seed=4893287):
    """
    Makes a realization of a naive SZ map

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * Number_of_SZ_Clusters: int, number of sources to add
        * Mean_Amplitude_of_SZ_Clusters: int, parameter for
            the exponential distribution
        * SZ_beta, SZ_Theta_core: floats, parameters for the beta function
    """
    SZMap = np.zeros((N,N))

    ## catalogue of SZ sources, X, Y, amplitude
    SZcat = np.zeros((3, Number_of_SZ_Clusters))
    state_initial = np.random.RandomState(seed)

    # make a distribution of point sources with varying amplitude
    i = 0
    while (i < Number_of_SZ_Clusters):
        pix_x = int(N * state_initial.rand())
        pix_y = int(N * state_initial.rand())
        pix_amplitude = state_initial.exponential(
            Mean_Amplitude_of_SZ_Clusters) * (-1.)
        SZcat[0,i] = pix_x
        SZcat[1,i] = pix_y
        SZcat[2,i] = pix_amplitude
        SZMap[pix_x, pix_y] += pix_amplitude
        i = i + 1

    if (do_plots):
        hist,bin_edges = np.histogram(
            SZMap, bins = 50, range=[SZMap.min(), -10])
        pl.semilogy(bin_edges[0:-1], hist)
        pl.xlabel('source amplitude [$\mu$K]')
        pl.ylabel('number or pixels')
        pl.show()

    # make a beta function
    beta = beta_function(N, pix_size, SZ_beta, SZ_Theta_core)

    # convovle the beta funciton with the point source amplitude to get the SZ map
    FT_beta = np.fft.fft2(np.fft.fftshift(beta))
    FT_SZMap = np.fft.fft2(np.fft.fftshift(SZMap))
    SZMap = np.fft.fftshift(np.real(np.fft.ifft2(FT_beta * FT_SZMap)))

    return SZMap, SZcat

def beta_function(N, pix_size, SZ_beta, SZ_Theta_core):
    """
    Make a beta function map

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * SZ_beta, SZ_Theta_core: floats, parameters for the beta function
    """
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) * pix_size
    X = np.outer(ones, inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)

    beta = (1 + (R / SZ_Theta_core)**2.)**((1 - 3. * SZ_beta) / 2.)

    return beta

def make_noise_map(
    N, pix_size,
    white_noise_level, atmospheric_noise_level,
    one_over_f_noise_level, seed=457349875):
    """
    Makes a realization of instrument noise, atmosphere and 1/f noise level
    set at 1 degrees

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * white_noise_level: float, level of white noise [uk.arcmin]
        * atmospheric_noise_level: float, level of atmospheric noise [uk.arcmin]
        * one_over_f_noise_level: float, level of 1/f noise [uk.arcmin]
    """
    ## make a white noise map
    state_initial = np.random.RandomState(seed)
    white_noise = state_initial.normal(
        0, 1, (N, N)) * white_noise_level / pix_size

    ## make an atmospheric noise map
    atmospheric_noise = 0.
    if (atmospheric_noise_level != 0):
        ones = np.ones(N)
        inds  = (np.arange(N) + .5 - N / 2.)
        X = np.outer(ones, inds)
        Y = np.transpose(X)

        ## angles relative to 1 degrees
        R = np.sqrt(X**2. + Y**2.) * pix_size /60.

        ## 0.01 is a regularization factor
        mag_k = 2 * np.pi/(R + .01)
        atmospheric_noise = np.fft.fft2(np.random.normal(0, 1, (N, N)))
        atmospheric_noise  = np.fft.ifft2(
            atmospheric_noise * np.fft.fftshift(
                mag_k**(5/3.))) * atmospheric_noise_level / pix_size

    ## make a 1/f map, along a single direction to illustrate striping
    oneoverf_noise = 0.
    if (one_over_f_noise_level != 0):
        ones = np.ones(N)
        inds  = (np.arange(N) + .5 - N/2.)

        ## angles relative to 1 degrees
        X = np.outer(ones, inds) * pix_size /60.
        ## 0.01 is a regularization factor
        kx = 2 * np.pi / (X + .01)
        oneoverf_noise = np.fft.fft2(np.random.normal(0, 1, (N, N)))
        oneoverf_noise = np.fft.ifft2(
            oneoverf_noise * np.fft.fftshift(kx)) * one_over_f_noise_level / pix_size

    ## return the noise map
    noise_map = np.real(white_noise + atmospheric_noise + oneoverf_noise)
    return noise_map

def cosine_window(N):
    """
    Makes a cosine window for apodizing to avoid edges effects in the 2d FFT

    Parameters
    -----------
        * N: int, number of pixels per row.
    """
    # make a 2d coordinate system
    ones = np.ones(N)

    ## eg runs from -pi/2 to pi/2
    inds  = (np.arange(N) + .5 - N / 2.) / N * np.pi
    X = np.outer(ones,inds)
    Y = np.transpose(X)

    # make a window map
    window_map = np.cos(X) * np.cos(Y)

    return window_map

def construct_noise_covar(
    N, pix_size, beam_size_fwhm,
    ell, DlTT,
    Number_of_Sources, Amplitude_of_Sources,
    Number_of_Sources_EX, Amplitude_of_Sources_EX,
    Number_of_SZ_Clusters, Mean_Amplitude_of_SZ_Clusters,
    SZ_beta, SZ_Theta_core,
    white_noise_level, atmospheric_noise_level, one_over_f_noise_level,
    n_iterations, seed=57467536):
    """
    Construct iteratively the 2d noise noise covariance in fourier space

    Parameters
    -----------
        * N: int, number of pixels per row
        * pix_size: float, pixel size
        * beam_size_fwhm: float, FWHM of the instrument
        * Number_of_Sources: int, number of faint sources to add
        * Amplitude_of_Sources: int, parameter for the Poisson distribution
        * Number_of_Sources_EX: int, number of bright sources to add
        * Amplitude_of_Sources_EX: int, parameter for the exponential distribution
        * Number_of_SZ_Clusters: int, number of sources to add
        * Mean_Amplitude_of_SZ_Clusters: int, parameter for
            the exponential distribution
        * SZ_beta, SZ_Theta_core: floats, parameters for the beta function
        * white_noise_level: float, level of white noise [uk.arcmin]
        * atmospheric_noise_level: float, level of atmospheric noise [uk.arcmin]
        * one_over_f_noise_level: float, level of 1/f noise [uk.arcmin]
        * n_iterations: int, number of iterations wanted to
            construct the covariance matrix


    """
    ## we will need a window function below to apodize the edge
    window = (cosine_window(N))

    FT_noise_covar = np.zeros((N, N))
    i = 0
    state_initial = np.random.RandomState(seed)
    seeds_it = state_initial.randint(1, 1e6, size=n_iterations)
    while (i < n_iterations):
        ## sumilate the astrophysical map
        CMB_T_temp = make_CMB_T_map(
            N, pix_size, ell, DlTT, seed=seeds_it[i])

        PSMap_temp = Poisson_source_component(
            N, pix_size, Number_of_Sources,
            Amplitude_of_Sources, seed=seeds_it[i])

        PSMap_temp += Exponential_source_component(
            N, pix_size, Number_of_Sources_EX,
            Amplitude_of_Sources_EX, seed=seeds_it[i])

        SZMap_temp, trash = SZ_source_component(
            N, pix_size, Number_of_SZ_Clusters,
            Mean_Amplitude_of_SZ_Clusters, SZ_beta,
            SZ_Theta_core, False, seed=seeds_it[i])

        CMB_T_temp  = CMB_T_temp + PSMap_temp + SZMap_temp

        ## fold in the instrument response
        CMB_T_convolved_temp = convolve_map_with_gaussian_beam(
            N, pix_size, beam_size_fwhm, CMB_T_temp)

        Noise_temp = make_noise_map(
            N, pix_size, white_noise_level,
            atmospheric_noise_level, one_over_f_noise_level)

        ## Fourier transform the map
        ## These are the two terms in the denominator
        temp =  np.fft.fft2(
            np.fft.fftshift( window * (CMB_T_convolved_temp + Noise_temp)))

        ## now average
        FT_noise_covar += np.real(
            np.conj(temp) * temp / (n_iterations * 1.0))

        ## note the progress
        sys.stdout.write(
            "\r matched filter noise realization, \
            iterations complete: %d of %d" % ((i + 1), n_iterations))
        sys.stdout.flush()

        ## iterate
        i = i + 1

    return FT_noise_covar

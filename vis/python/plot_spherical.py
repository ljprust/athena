#! /usr/bin/env python

"""
Script for plotting vertical (r,theta) or midplane (r,phi) slices of data in
spherical coordinates.

Run "plot_spherical.py -h" to see description of inputs.

See documentation on athena_read.athdf() for important notes about reading files
with mesh refinement.

Users are encouraged to make their own versions of this script for improved
results by adjusting figure size, spacings, tick locations, axes labels, etc.
The script must also be modified to plot any functions of the quantities in the
file, including combinations of multiple quantities.

Requires scipy if making a stream plot.
"""

# Python standard modules
import argparse
import warnings

# Other Python modules
import numpy as np

# Athena++ modules
import athena_read

# for parsing text
import re

plotCircle    = False      # plots a circle centered on the origin
circleRadius  = 7.0e9
nfiles        = 1          # number of files to analyze
first         = 100          # number of first file
interval      = 1          # interval between file numbers
fileprefix    = 'wt.out1.' # file name before number
filesuffix    = '.athdf'   # file name after number
paramfilename = None # 'athinput.windtunnel_wedge'

# pattern matching the line "GM = ?" in parameter file (can be in scientific notation)
GM_pattern = r"^GM\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
GM_value = None

gamma_pattern = r"^gamma\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
gamma = None

print('Parsing parameter file...')
if (paramfilename is not None) :
    with open(paramfilename, "r") as file:
        for line in file:
            match = re.search(gamma_pattern, line)
            if match:
                # extract the string representation and convert to a float
                gamma = float(match.group(1))
                break # stop checking the file once found
if (paramfilename is not None) :
    # open and process the parameter file
    with open(paramfilename, "r") as file:
        for line in file:
            match = re.search(GM_pattern, line)
            if match:
                # extract the string representation and convert to a float
                GM_value = float(match.group(1))
                break # stop checking the file once found

# verify and use the float value
if GM_value is not None:
    print(f"    Detected GM = {GM_value}")
else:
    GM_value = 0.0
    print("    Did not read GM from parameter file")

if gamma is not None:
    print(f"    Detected gamma = {gamma}")
else:
    gamma = 0.0
    print("    Did not read gamma from parameter file")

# create list of filenames
filename = []
for i in range(0, nfiles) :
    num = 100000 + first + i * interval
    numstr = str(num)
    cut = numstr[1:7]
    filename.append(cut)

# Main function
def main(myj,**kwargs):

    if(nfiles>0):       
        kwargs['data_file'] = fileprefix+filename[myj]+filesuffix
    print('starting file',kwargs['data_file'])

    # Load function for transforming coordinates
    if kwargs['stream'] is not None:
        from scipy.ndimage import map_coordinates

    # Load Python plotting modules
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    # Determine refinement level to use
    if kwargs['level'] is not None:
        level = kwargs['level']
    else:
        level = None

    # Determine if vector quantities should be read
    quantities = [kwargs['quantity']]
    if kwargs['stream'] is not None:
        quantities.append(kwargs['stream'] + '1')
        if kwargs['midplane']:
            quantities.append(kwargs['stream'] + '3')
        else:
            quantities.append(kwargs['stream'] + '2')

    # Define grid compression in theta-direction
    h = kwargs['theta_compression'] if kwargs['theta_compression'] is not None else 1.0

    def theta_func(xmin, xmax, _, nf):
        x2_vals = np.linspace(xmin, xmax, nf)
        theta_vals = x2_vals + (1.0 - h) / 2.0 * np.sin(2.0 * x2_vals)
        return theta_vals

    # Read data
    if kwargs['theta_compression'] is not None:
        if quantities[0] == 'Levels':
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities[1:],
                                     level=level, return_levels=True,
                                     face_func_2=theta_func)
        else:
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level, face_func_2=theta_func)
    else:
        if quantities[0] == 'Levels':
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities[1:],
                                     level=level, return_levels=True)
        elif quantities[0] == 'entropy' :
            cbar_label = r'$P/\rho^\gamma$'
            saveasprefix = 'entropy'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['entropy'] = data['press']*np.power(data['rho'],-gamma)
        elif quantities[0] == 'enthalpy' :
            cbar_label = r'$h$'
            saveasprefix = 'enthalpy'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['enthalpy'] = gamma/(gamma-1.0)*data['press']/data['rho']
        elif quantities[0] == 'r0' :
            cbar_label = r'$\chi$'
            saveasprefix = 'r0_'
            data = athena_read.athdf(kwargs['data_file'], quantities=['r0'],
                                     level=level)
        elif quantities[0] == 'rhor0' :
            cbar_label = r'$\chi\rho$'
            saveasprefix = 'rhor0_'
            data = athena_read.athdf(kwargs['data_file'], level=level)
            data['rhor0'] = data['rho']*data['r0']
        elif quantities[0] == 'bernoulli' :
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['bernoulli'] = 0.5*data['rho']*(data['vel1']*data['vel1']+data['vel2']*data['vel2']+data['vel3']*data['vel3']) + data['press']
        elif quantities[0] == 'energy' :
            cbar_label = 'Specific Energy'
            saveasprefix = 'energy'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['energy'] = 0.5*(data['vel1']*data['vel1']+data['vel2']*data['vel2']+data['vel3']*data['vel3']) + data['press']/(gamma-1.0)/data['rho']
        elif quantities[0] == 'kinetic' :
            cbar_label = 'Kinetic Energy'
            saveasprefix = 'kinetic'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['kinetic'] = 0.5*(data['vel1']*data['vel1']+data['vel2']*data['vel2']+data['vel3']*data['vel3'])
        elif quantities[0] == 'totalenthalpy' :
            cbar_label = r'$h_{\rm tot}$'
            saveasprefix = 'htot'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            data['totalenthalpy'] = 0.5*(data['vel1']*data['vel1']+data['vel2']*data['vel2']+data['vel3']*data['vel3']) + data['press']*gamma/(gamma-1.0)/data['rho']
        elif quantities[0] == 'mach' :
            cbar_label = r'$\mathcal{M}$'
            saveasprefix = 'mach'
            data      = athena_read.athdf(kwargs['data_file'],
                                     level=level)
            cs2 = gamma*data['press']/data['rho']
            v2 = data['vel1']*data['vel1']+data['vel2']*data['vel2']+data['vel3']*data['vel3']
            data['mach'] = np.sqrt(v2/cs2)
            '''
            elif kwargs['vorticity'] :
                data = athena_read.athdf(kwargs['data_file'], quantities=['vel1','vel2','vel3'], level=level)
                coordinates = data['Coordinates'].decode('ascii', 'replace')
                r = data['x1v']
                theta = data['x2v']
                phi = data['x3v']
                r_face = data['x1f']
                nx1 = len(r)
                nx2 = len(theta)
                nx3 = len(phi)           
                Pgradfront[l] = (pressfront[l+1]-pressfront[l])/(r[l+1]-r[l])
                    Pgradrhofront[l] = 0.5*(rhofront[l]+rhofront[l+1])
                    Pgradback[l] = (pressback[l+1]-pressback[l])/(r[l+1]-r[l])
                    Pgradrhoback[l] = 0.5*(rhoback[l]+rhoback[l+1])
                    Pgradrad[l] = r_face[l+1]
            '''
        elif quantities[0] == 'rho' :
            cbar_label = r'$\rho$'
            saveasprefix = 'rho'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'press' :
            cbar_label = r'$P$'
            saveasprefix = 'pres'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel1' :
            cbar_label = r'$v_r$'
            saveasprefix = 'vr'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel2' :
            cbar_label = r'$v_\theta$'
            saveasprefix = 'vtheta'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel3' :
            cbar_label = r'$v_\phi$'
            saveasprefix = 'vphi'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel_xyz1' :
            cbar_label = r'$v_x$'
            saveasprefix = 'vx'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel_xyz2' :
            cbar_label = r'$v_y$'
            saveasprefix = 'vy'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        elif quantities[0] == 'vel_xyz3' :
            cbar_label = r'$v_z$'
            saveasprefix = 'vz'
            data = athena_read.athdf(kwargs['data_file'], quantities=quantities,
                                     level=level)
        else :
            printf('QUANTITY NOT DETECTED')
        time = data['Time']

    # Extract basic coordinate information
    coordinates = data['Coordinates'].decode('ascii', 'replace')
    #r = data['x1v'] # this gives NaNs for some reason
    theta = data['x2v']
    phi = data['x3v']
    r_face = data['x1f']
    theta_face = data['x2f']
    phi_face = data['x3f']
    nx1 = len(r_face-1)
    nx2 = len(theta)
    nx3 = len(phi)
    r = 0.5*(r_face[1:nx1]+r_face[0:nx1-1])

    # Set radial extent
    if kwargs['r_max'] is not None:
        r_max = kwargs['r_max']
    else:
        r_max = r_face[-1]

    # Account for logarithmic radial coordinate
    if kwargs['logr']:
        r = np.log10(r)
        r_face = np.log10(r_face)
        r_max = np.log10(r_max)

    # Create scalar grid
    if kwargs['midplane']:
        r_grid, phi_grid = np.meshgrid(r_face, phi_face)
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)
    else:
        theta_face_extended = np.concatenate((theta_face, 2.0*np.pi - theta_face[-2::-1]))
        r_grid, theta_grid = np.meshgrid(r_face, theta_face_extended)
        x_grid = r_grid * np.sin(theta_grid)
        y_grid = r_grid * np.cos(theta_grid)

    # Create streamline grid
    if kwargs['stream'] is not None:
        x_stream = kwargs['xoffset']*kwargs['lscale']+np.linspace(-r_max, r_max, kwargs['stream_samples'])
        if kwargs['midplane']:
            y_stream = np.linspace(-r_max, r_max, kwargs['stream_samples'])
            x_grid_stream, y_grid_stream = np.meshgrid(x_stream, y_stream)
            r_grid_stream_coord = (x_grid_stream.T**2 + y_grid_stream.T**2) ** 0.5
            phi_grid_stream_coord = np.pi + np.arctan2(-y_grid_stream.T, -x_grid_stream.T)
            phi_grid_stream_pix = ((phi_grid_stream_coord + phi[0])
                                   / (2.0*np.pi + 2.0 * phi[0])) * (nx3 + 1)
        else:
            z_stream = np.linspace(-r_max, r_max, kwargs['stream_samples'])
            x_grid_stream, z_grid_stream = np.meshgrid(x_stream, z_stream)
            r_grid_stream_coord = (x_grid_stream.T**2 + z_grid_stream.T**2) ** 0.5
            theta_grid_stream_coord = np.pi - \
                np.arctan2(x_grid_stream.T, -z_grid_stream.T)
            if kwargs['theta_compression'] is None:
                theta_grid_stream_pix = ((theta_grid_stream_coord + theta[0])
                                         / (2.0*np.pi + 2.0 * theta[0])) * (2 * nx2 + 1)
            else:
                theta_grid_stream_pix = np.empty_like(theta_grid_stream_coord)
                theta_extended = np.concatenate((-theta[0:1], theta,
                                                 2.0*np.pi - theta[::-1],
                                                 2.0*np.pi + theta[0:1]))
                for (i, j), theta_val in np.ndenumerate(theta_grid_stream_coord):
                    index = sum(theta_extended[1:-1] < theta_val) - 1
                    if index < 0:
                        theta_grid_stream_pix[i, j] = -1
                    elif index < 2 * nx2 - 1:
                        theta_grid_stream_pix[i, j] = (
                            index + ((theta_val - theta_extended[index])
                                     / (theta_extended[index+1] - theta_extended[index])))
                    else:
                        theta_grid_stream_pix[i, j] = 2 * nx2 + 2
        r_grid_stream_pix = np.empty_like(r_grid_stream_coord)
        for (i, j), r_val in np.ndenumerate(r_grid_stream_coord):
            index = sum(r < r_val) - 1
            if index < 0:
                r_grid_stream_pix[i, j] = -1
            elif index < nx1 - 1:
                r_grid_stream_pix[i, j] = index + \
                    (r_val - r[index]) / (r[index + 1] - r[index])
            else:
                r_grid_stream_pix[i, j] = nx1

    # Perform slicing/averaging of scalar data
    if kwargs['midplane']:
        if nx2 % 2 == 0:
            vals = np.mean(data[kwargs['quantity']][:, nx2//2-1:nx2//2+1, :], axis=1)
        else:
            vals = data[kwargs['quantity']][:, nx2//2, :]
        if kwargs['average']:
            vals = np.repeat(np.mean(vals, axis=0, keepdims=True), nx3, axis=0)
    else:
        if kwargs['average']:
            vals_right = np.mean(data[kwargs['quantity']], axis=0)
            vals_left = vals_right
        else:
            vals_right = 0.5 * (data[kwargs['quantity']]
                                [-1, :, :] + data[kwargs['quantity']][0, :, :])
            vals_left = 0.5 * (data[kwargs['quantity']][(nx3//2)-1, :, :]
                               + data[kwargs['quantity']][nx3//2, :, :])
        '''
        else:
            vals_right = 0.5 * (data[kwargs['quantity']][nx3/4-1, :, :]
                               + data[kwargs['quantity']][nx3/4, :, :])
            vals_left = 0.5 * (data[kwargs['quantity']][nx3*3/4-1, :, :]
                               + data[kwargs['quantity']][nx3*3/4, :, :])
        '''

    # Join scalar data through boundaries
    if not kwargs['midplane']:
        vals = np.vstack((vals_right, vals_left[::-1, :]))

    # Perform slicing/averaging of vector data
    if kwargs['stream'] is not None:
        if kwargs['midplane']:
            if nx2 % 2 == 0:
                vals_r = np.mean(data[kwargs['stream'] + '1']
                                 [:, nx2//2-1:nx2//2+1, :], axis=1).T
                vals_phi = np.mean(data[kwargs['stream'] + '3']
                                   [:, nx2//2-1:nx2//2+1, :], axis=1).T
            else:
                vals_r = data[kwargs['stream'] + '1'][:, nx2//2, :].T
                vals_phi = data[kwargs['stream'] + '3'][:, nx2//2, :].T
            if kwargs['stream_average']:
                vals_r = np.tile(np.reshape(np.mean(vals_r, axis=1), (nx1, 1)), nx3)
                vals_phi = np.tile(np.reshape(np.mean(vals_phi, axis=1), (nx1, 1)), nx3)
        else:
            if kwargs['stream_average']:
                vals_r_right = np.mean(data[kwargs['stream'] + '1'], axis=0).T
                vals_r_left = vals_r_right
                vals_theta_right = np.mean(data[kwargs['stream'] + '2'], axis=0).T
                vals_theta_left = -vals_theta_right
            else:
                vals_r_right = data[kwargs['stream'] + '1'][0, :, :].T
                vals_r_left = data[kwargs['stream'] + '1'][nx3//2, :, :].T
                vals_theta_right = data[kwargs['stream'] + '2'][0, :, :].T
                vals_theta_left = -data[kwargs['stream'] + '2'][nx3//2, :, :].T

    # Join vector data through boundaries
    if kwargs['stream'] is not None:
        if kwargs['midplane']:
            vals_r = np.hstack((vals_r[:, -1:], vals_r, vals_r[:, :1]))
            vals_r = map_coordinates(vals_r, (r_grid_stream_pix, phi_grid_stream_pix),
                                     order=1, cval=np.nan)
            vals_phi = np.hstack((vals_phi[:, -1:], vals_phi, vals_phi[:, :1]))
            vals_phi = map_coordinates(vals_phi, (r_grid_stream_pix, phi_grid_stream_pix),
                                       order=1, cval=np.nan)
        else:
            vals_r = np.hstack((vals_r_left[:, :1], vals_r_right, vals_r_left[:, ::-1],
                                vals_r_right[:, :1]))
            vals_r = map_coordinates(vals_r, (r_grid_stream_pix, theta_grid_stream_pix),
                                     order=1, cval=np.nan)
            vals_theta = np.hstack((vals_theta_left[:, :1], vals_theta_right,
                                    vals_theta_left[:, ::-1], vals_theta_right[:, :1]))
            vals_theta = map_coordinates(vals_theta,
                                         (r_grid_stream_pix, theta_grid_stream_pix),
                                         order=1, cval=np.nan)

    # Transform vector data to Cartesian components
    if kwargs['stream'] is not None:
        if kwargs['logr']:
            r_vals = 10.0**r_grid_stream_coord
            logr_vals = r_grid_stream_coord
        else:
            r_vals = r_grid_stream_coord
        if kwargs['midplane']:
            sin_phi = np.sin(phi_grid_stream_coord)
            cos_phi = np.cos(phi_grid_stream_coord)
            if kwargs['logr']:
                dx_dr = 1.0 / (np.log(10.0) * r_vals) * cos_phi
                dy_dr = 1.0 / (np.log(10.0) * r_vals) * sin_phi
                dx_dphi = -logr_vals * sin_phi
                dy_dphi = logr_vals * cos_phi
            else:
                dx_dr = cos_phi
                dy_dr = sin_phi
                dx_dphi = -r_vals * sin_phi
                dy_dphi = r_vals * cos_phi
            if not (coordinates == 'schwarzschild' or coordinates == 'kerr-schild'):
                dx_dphi /= r_vals
                dy_dphi /= r_vals
            vals_x = dx_dr * vals_r + dx_dphi * vals_phi
            vals_y = dy_dr * vals_r + dy_dphi * vals_phi
        else:
            sin_theta = np.sin(theta_grid_stream_coord)
            cos_theta = np.cos(theta_grid_stream_coord)
            if kwargs['logr']:
                dx_dr = 1.0 / (np.log(10.0) * r_vals) * sin_theta
                dz_dr = 1.0 / (np.log(10.0) * r_vals) * cos_theta
                dx_dtheta = logr_vals * cos_theta
                dz_dtheta = -logr_vals * sin_theta
            else:
                dx_dr = sin_theta
                dz_dr = cos_theta
                dx_dtheta = r_vals * cos_theta
                dz_dtheta = -r_vals * sin_theta
            if not (coordinates == 'schwarzschild' or coordinates == 'kerr-schild'):
                dx_dtheta /= r_vals
                dz_dtheta /= r_vals
            vals_x = dx_dr * vals_r + dx_dtheta * vals_theta
            vals_z = dz_dr * vals_r + dz_dtheta * vals_theta

    # Determine colormapping properties
    cmap = plt.get_cmap(kwargs['colormap'])
    vmin = kwargs['vmin']
    vmax = kwargs['vmax']
    if kwargs['logc']:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

    # Make plot
    plt.figure()
    im = plt.pcolormesh(x_grid/kwargs['lscale'], y_grid/kwargs['lscale'], vals, cmap=cmap, norm=norm)
    #cont = plt.contour(y_grid[0:-1,0:-1], x_grid[0:-1,0:-1], vals, 1, colors='k',origin='lower')
    plt.gca().set_aspect('equal')
    plt.gca().set_ylim(-r_max/kwargs['lscale']+kwargs['xoffset'], r_max/kwargs['lscale']+kwargs['xoffset'])
    plt.gca().set_xlim(-r_max/kwargs['lscale'], r_max/kwargs['lscale'])
    plt.xticks(fontsize=14) # 15)
    plt.yticks(fontsize=14) # 15)
    if plotCircle :
        circle = plt.Circle((0.0,0.0), circleRadius, fc='None', ec='r', lw=1.0, ls='--')
        plt.gca().add_patch(circle)
    if kwargs['stream'] is not None:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                'invalid value encountered in greater_equal',
                RuntimeWarning,
                'numpy')
            if kwargs['midplane']:
                plt.streamplot(x_stream/kwargs['lscale'], y_stream/kwargs['lscale'], vals_x.T, vals_y.T,
                               density=kwargs['stream_density'], color='k')
            else:
                plt.streamplot(x_stream/kwargs['lscale'], z_stream/kwargs['lscale'], vals_x.T, vals_z.T,
                               density=kwargs['stream_density'], color='k')
    plt.gca().set_aspect('equal')
    if kwargs['logr']:
        if kwargs['midplane']:
            plt.xlabel(r'$\log_{10}(r)\ x / r$')
            plt.ylabel(r'$\log_{10}(r)\ y / r$')
        else:
            plt.xlabel(r'$\log_{10}(r)\ x / r$')
            plt.ylabel(r'$\log_{10}(r)\ z / r$')
    else:
        if kwargs['midplane']:
            plt.xlabel(r'$x$')
            plt.ylabel(r'$y$')
        else:
            plt.xlabel(r'$x$')
            plt.ylabel(r'$z$')
    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label, fontsize=16)
    #myboundaries = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
    #cbar = plt.colorbar(im, shrink=0.887, anchor=(0.0,0.78), boundaries=myboundaries) # 0.8
    #ax = plt.axes()
    #cbar = plt.colorbar(im, boundaries=myboundaries)
    #cbar.ax.tick_params(labelsize=14)

    #plt.annotate( r'$t = $' + str(time)[0:4] + r' $R_{A}/c_{s}$', xy=(0.02,0.02), xytext=(0.02,0.02), xycoords='axes fraction', color='white', size=16 ) # 16
    #plt.annotate( r'$\mathcal{M}_{\infty}=0.6$, $\gamma=5/3$', xy=(0.02,0.02), xytext=(0.02,0.02), xycoords='axes fraction', color='red', size=20 )
    #ax.text(0.02, 0.02, r'$\mathcal{M}_{\infty}=0.6$, $\gamma=5/3$', color='black', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='black'), size=16)
    #plt.tight_layout()

    if(nfiles>0):
        kwargs['output_file'] = saveasprefix + filename[myj]
    if kwargs['output_file'] == 'show':
        plt.show()
    else:
        plt.savefig(kwargs['output_file'], bbox_inches='tight')
        print('saved figure',kwargs['output_file'])


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('data_file',
    #                    help='name of input file, possibly including path')
    parser.add_argument('quantity',
                        help='name of quantity to be plotted')
    #parser.add_argument('output_file',
    #                    help=('name of output to be (over)written, possibly including '
    #                          'path; use "show" to show interactive plot instead'))
    parser.add_argument('-m',
                        '--midplane',
                        action='store_true',
                        help=('flag indicating plot should be midplane (r,phi) rather '
                              'than (r,theta)'))
    parser.add_argument('-a', '--average',
                        action='store_true',
                        help='flag indicating phi-averaging should be done')
    parser.add_argument('-l',
                        '--level',
                        type=int,
                        default=None,
                        help=('refinement level to be used in plotting (default: max '
                              'level in file)'))
    parser.add_argument('-r', '--r_max',
                        type=float,
                        default=None,
                        help='maximum radial extent of plot')
    parser.add_argument('--logr',
                        action='store_true',
                        help='flag indicating data should be plotted logarithmically in '
                             'radius')
    parser.add_argument('-c',
                        '--colormap',
                        default=None,
                        help=('name of Matplotlib colormap to use instead of default'))
    parser.add_argument('--vmin',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap minimum; use '
                              '--vmin=<val> if <val> has negative sign'))
    parser.add_argument('--vmax',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap maximum; use '
                              '--vmax=<val> if <val> has negative sign'))
    parser.add_argument('--logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')
    parser.add_argument('-s', '--stream',
                        default=None,
                        help='name of vector quantity to use to make stream plot')
    parser.add_argument('--stream_average',
                        action='store_true',
                        help='flag indicating phi-averaging on stream plot data')
    parser.add_argument('--stream_density',
                        type=float,
                        default=1.0,
                        help='density of stream lines')
    parser.add_argument('--stream_samples',
                        type=int,
                        default=100,
                        help='linear size of stream line sampling grid')
    parser.add_argument('--theta_compression',
                        type=float,
                        default=None,
                        help=('compression parameter h in '
                              'theta = pi*x_2 + (1-h)/2 * sin(2*pi*x_2)'))
    parser.add_argument('--xoffset',
                        type=float,
                        default=0.0,
                        help='x offset of plot center')
    parser.add_argument('--lscale',
                        type=float,
                        default=1.0,
                        help=('length scale'))
    args = parser.parse_args()
    jindex=0
    for name in filename :
        main(jindex,**vars(args))
        jindex = jindex+1 

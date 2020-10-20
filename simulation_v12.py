#!/usr/bin/env python3
# coding=utf-8

import sys
from datetime import datetime
import argparse
import pandas as pd
import numpy as np
import math
from matplotlib import cm
from matplotlib import collections as mc
import copy
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import PatchCollection
import imageio as io
import os
import matplotlib.colors as clr
matplotlib.rcParams.update({'font.size': 6})


def get_parser():
    # build commandline parser
    parser = argparse.ArgumentParser(
        description="MYCEMULATOR – Simulator of mycelial growth\n"
                    "As you can see, there are many parameters to tune in order to make the simulation fit your specific purpose.\n"
                    "Default values are all based on the growth of ATCC 1015.\n"
                    "Enjoy!")
    parser.add_argument("-ani","--animation", type=str, dest="ani", default=None, help="Name of folder of time-dependent png images and gif animation."
                                                                                 "An mp4 can easily be created from the images by using ffmpeg")
    parser.add_argument("-pdf", type=str, dest="pdf", default=None, help="Name of pdf of time-dependent subplots.")
    parser.add_argument("-img", "--img_ana", type=str, dest="img_ana",
                        help="Use parameters from csv-output of image analysis tool")
    parser.add_argument("--hours", type=float, dest="hours", default=12,
                        help="Number of hours to simulate exponential hyphal growth for. Default: 12")
    parser.add_argument("-tstep", type=float, dest="tstep", default=1/60,
                        help="Number/part of hours per simulation round. Default: 1/60, i.e. 1 minute per round")
    parser.add_argument("-field", "--fieldtype", type=str, dest="fieldtype", default="g", choices=['g', 'u'],
                        help="Substrate field type. Choose whether to use a uniform field (u) or a radial gradient field (g).")
    parser.add_argument("-source", action='store_true', dest="source", default=False,
                        help="Use to have an infinite substrate source at the edge of the simulation area.")
    parser.add_argument("-mu_max", type=float, dest="mu_max", default=0.275733333333333,
                        help="Mu_max. Maximum growth rate for the given strain in exponential phase. Default is biolector measurement of ATCC 1015")
    parser.add_argument("-q", type=float, dest="q", default=0.016925731922538392,
                        help="Branching frequency per hyphal element per hour (if not taken from image analysis output). Default is taken from image analysis of ATCC 1015")
    parser.add_argument("-lat_sub_min", type=float, dest="lat_sub_min", default=10,
                        help="Minimum substrate concentration for lateral branching to occur")
    parser.add_argument("-ap_sub_min", type=float, dest="ap_sub_min", default=10,
                        help="Minimum substrate concentration for apical branching to occur")
    parser.add_argument("-ext_sub_min", type=float, dest="ext_sub_min", default=5,
                        help="Minimum substrate concentration for tip extension to occur")
    parser.add_argument("-S0", type=float, dest="S0", default=100,
                        help="Initial uniform substrate concentration. Also boundary condition if '-source' is applied.")
    parser.add_argument("-D", type=float, dest="D", default=0.99,
                        help="Diffusion constant. Default: 0.99")
    parser.add_argument("-S_tip", type=float, dest="S_tip", default=10,
                        help="Substrate consumption of a hyphal tip")
    parser.add_argument("-S_nontip", type=float, dest="S_nontip", default=1,
                        help="Substrate consumption of a non-tip hyphal element")
    return parser


def get_args():
    parser = get_parser()
    args = parser.parse_args()
    return args


def img_ana_parameters(csv_file):
    return pd.read_csv(csv_file, header=0)


def tip_extension_monod(mu_max, hyph_row, hyphal_elements, St, Ks, ext_sub_min, curvature_gamma_params, r, time):
    """source: Lejeune et al 1995, Morphology of Trichoderma reesei QM 9414 in Submerged Cultures"""
    sep_distance = center_distance(hyph_row)
    if sep_distance >= 20:
        return hyphal_elements

    extension = mu_max * St[int(r/20*sep_distance)] / (St[int(r/20*sep_distance)] + Ks)
    if St[int(r/20*sep_distance)] > ext_sub_min and extension > np.random.uniform(0, 1):
        # get coordinates of old tip
        x0_old, y0_old, x1_old, y1_old = midpoints_to_endpoints(hyph_row)
        if 0 < hyph_row.angle < 180:
            x_oldtip, y_oldtip = max(((x0_old, y0_old), (x1_old, y1_old)),  key=lambda x: x[1])
        elif 180 < hyph_row.angle < 360:
            x_oldtip, y_oldtip = min(((x0_old, y0_old), (x1_old, y1_old)),  key=lambda x: x[1])
        elif hyph_row.angle == 0:
            x_oldtip, y_oldtip = max(((x0_old, y0_old), (x1_old, y1_old)), key=lambda x: x[0])
        elif hyph_row.angle == 180:
            x_oldtip, y_oldtip = min(((x0_old, y0_old), (x1_old, y1_old)), key=lambda x: x[0])

        new_coord, center_dist, angles = list(), list(), list()
        for _ in range(10):
            new_angle = hyph_row.angle + np.random.choice((-1, 1)) * round(5*np.random.gamma(*curvature_gamma_params))
            # check that new angle is >=0 and < 360
            if not 0 <= new_angle < 360:
                new_angle = new_angle % 360

            # find out if x and y positions for new branch are smaller/larger than before
            if 0 < new_angle < 180:
                y_dir = 1
            else:
                y_dir = -1
            if 0 < new_angle < 90 or 270 < new_angle < 360:
                x_dir = 1
            else:
                x_dir = -1

            # make sure you get right angle of triangle
            if 0 < new_angle < 90: use_angle = new_angle
            elif 90 < new_angle < 180:
                use_angle = 180 - new_angle
            elif 180 < new_angle < 270:
                use_angle = new_angle - 180
            elif 270 < new_angle < 360:
                use_angle = 360 - new_angle

            if new_angle in (0, 90, 180, 270):
                if new_angle in (0, 180):
                    delta_x = x_dir * 0.5
                    delta_y = 0
                elif new_angle in (90, 270):
                    delta_x = 0
                    delta_y = y_dir * 0.5

            else:
                # change in y direction based on sinus relations
                delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2 * np.sin(np.deg2rad(90))))
                # change in x direction based on normal pythagoras
                delta_x = x_dir * np.sqrt(0.5 ** 2 - delta_y ** 2)

            x_mid, y_mid = x_oldtip + delta_x, y_oldtip + delta_y

            new_coord.append((x_mid, y_mid))
            dist = np.sqrt(x_mid ** 2 + y_mid ** 2)
            if dist > 20: dist = 20
            center_dist.append(dist)
            angles.append(new_angle)

            # if substrate concentration is more than 1% as high for one proposed angle than the first,
            # choose this angle, else choose first
            if St[int(r/20*max(center_dist))] > 1.01 * St[int(r/20*center_dist[0])]:
                print('tropism!')
                new_x_mid, new_y_mid = new_coord[np.argmax(center_dist)]
                new_angle = angles[np.argmax(center_dist)] % 360

            else:
                new_x_mid, new_y_mid = new_coord[0]
                new_angle = angles[0] % 360

        hyphal_elements.loc[len(hyphal_elements)] = {'x_mid':new_x_mid, 'y_mid':new_y_mid, 'angle': new_angle, 'tip': True, 'time':time}

        # change tip status of hyph_row, as this hyphal element has been extended and is no longer a tip
        hyphal_elements.at[hyph_row.name, 'tip'] = False
    return hyphal_elements


def hyphal_length(hyphal_elements):
  return math.sqrt((hyphal_elements.x-hyphal_elements.x0)**2 + (hyphal_elements.y-hyphal_elements.y0)**2)


def substrate_per_distance(n_tip, n_nontip, St, S0=200, r=200, h=1, D=0.1, S_tip=5, S_nontip=5, fieldtype = "g", source = False):
    """
    :param n_tip: number of tip hyphal_elements
    :param n_nontip: number of non-tip hyphal_elements
    :param St: St list of substrate values for previous simulation step
    :param S0: initial substrate concentration - also boundary condition
    :param r: radius of substrate field
    :param h: radial step size
    :param D: diffusion constant
    :param S_tip: amount of substrate used per tip hyphal element
    :param S_nontip: amount of substrate used per non-tip hyphal element
	:param fieldtype: choose gradient or uniform field
	:param source: choose whether to have an infinite source at the edge of the field or not
    :return: St: list of increasing substrate concentrations as a function of distance to center
    """

    # uniform substrate field
    if fieldtype == "u":
        if source is False:
            if St[0] == 0:
                dSdt = [0 for i in range(r+1)]
            else:
                dSdt = [- (S_tip * n_tip + S_nontip * n_nontip)/(St[0]*r) for i in range(r+1)]
        else:
            dSdt = [- (S_tip * n_tip + S_nontip * n_nontip)/(St[0]*r) for i in range(r)]
            dSdt.append(0)

    # gradial substrate field
    elif fieldtype == "g":
        dSdt = [0 for i in range(r + 1)]  # list to hold changes in substrate levels for this timestep
        # continuous equation is dS/dt = D*d2S/dr2 (minus use at centre) - using difference estimates for second order derivative
		# forward difference for center
        dSdt[0] = min(0, D * (St[0 + 2 * h] - 2 * St[1] + St[0]) / h ** 2 - (S_tip * n_tip + S_nontip * n_nontip))
        for rad in range(1, r):
            # centre difference for most of the field - also discarding positive differences which happen due to numerical solution
            dSdt[rad] = min(0, D * (St[rad + 1] - 2 * St[rad] + St[rad - 1]) / h ** 2)
        # centre difference but with S(rad>r)=S0
        if source is False:
            dSdt[r] = min(0, D * (-St[r] + St[r - 1]) / h ** 2)

    # add change in substrate level to previous level
    St = [max(0, a + b) for a, b in zip(St, dSdt)]

    # make sure substrate levels do not exceed initial level
    if max(St) > S0:
        St = [min(st, S0) for st in St]

    return St


def center_distance(hyph_row):
    """Compute distance to center for a given hyphal element"""
    return np.sqrt(hyph_row.x_mid**2+hyph_row.y_mid**2)


def branching(hyph_row, St, min_substrate, p_branching, hyphal_elements, branch_angle_beta_params, r, time):
    """
    Add laterally or apically branched hyphal_elements (by a certain probability) to hyphal_elements if enough substrate
    :param hyph_row: row from hyphal_elements df
    :param St: substrate concentration per distance for this timestep
    :param min_substrate: minimum substrate for branching to occur
    :param p_branching: probability of branching event, given enough substrate
    :param hyphal_elements: hyphal_elements df to add evt. branched septum to
    :return:
    """
    # if given hyphal element hits the edge of the field, no branching event
    sep_distance = center_distance(hyph_row)
    if sep_distance >= 20:
        return hyphal_elementss

    # branching occurs at a given probability if enough substrate
    if St[int(r/20*sep_distance)] > min_substrate and p_branching > np.random.uniform(0, 1):
        new_angle = hyph_row.angle + np.random.choice((-1,1)) * round(np.random.beta(*branch_angle_beta_params)*90)      # 90 is the largest branching angle
        # check that new angle does not exceed 360
        if not 0 <= new_angle < 360: new_angle = new_angle % 360

        # find out if x and y positions for new branch are smaller/larger than before
        if 0 < new_angle < 180: y_dir = 1
        else: y_dir = -1
        if 0 < new_angle < 90 or 270 < new_angle < 360: x_dir = 1
        else: x_dir = -1

        # make sure you get right angle of triangle
        if 0 < new_angle < 90: use_angle = new_angle
        if 90 < new_angle < 180: use_angle = 180 - new_angle
        elif 180 < new_angle < 270: use_angle = new_angle - 180
        elif 270 < new_angle < 360: use_angle = 360 - new_angle

        if new_angle in (0, 90, 180, 270):
            if new_angle in (0,180):
                delta_x = x_dir*0.5
                delta_y = 0
            elif new_angle in (90,270):
                delta_x = 0
                delta_y = y_dir*0.5

        else:
            # change in y direction based on sinus relations
            delta_y = y_dir * abs(np.sin(np.deg2rad(use_angle)) / (2*np.sin(np.deg2rad(90))))
            # change in x direction based on normal pythagoras
            delta_x = x_dir * np.sqrt(0.5**2-delta_y**2)

        new_x_mid, new_y_mid = hyph_row.x_mid + delta_x, hyph_row.y_mid + delta_y

        hyphal_elements.loc[len(hyphal_elements)] = {'x_mid':new_x_mid, 'y_mid':new_y_mid, 'angle': new_angle, 'tip': True, 'time':time}

    return hyphal_elements


def midpoints_to_endpoints(hyph_row):
    """
    Use x_mid, y_mid and angle to find end coordinates of hyphal_elements
    """

    # return for easy angles
    if hyph_row.angle in (0, 90, 180, 270):
        if hyph_row.angle in (0,180):
            return hyph_row.x_mid-0.5, hyph_row.y_mid, hyph_row.x_mid+0.5, hyph_row.y_mid
        elif hyph_row.angle in (90,270):
            return hyph_row.x_mid, hyph_row.y_mid - 0.5, hyph_row.x_mid, hyph_row.y_mid + 0.5

    # get right angle of triangle to use for sinus relations
    if hyph_row.angle < 90: use_angle = hyph_row.angle
    if 90 < hyph_row.angle < 180: use_angle = 180 - hyph_row.angle
    if 180 < hyph_row.angle < 270: use_angle = hyph_row.angle - 180
    if 270 < hyph_row.angle < 360: use_angle = 360 - hyph_row.angle

    # change in y direction based on sinus relations
    delta_y = abs(np.sin(np.deg2rad(use_angle)) / (2 * np.sin(np.deg2rad(90))))
    # change in x direction based on normal pythagoras
    delta_x = np.sqrt(abs(0.5 ** 2 - delta_y ** 2))

    # check how to add delta_x and delta_y
    if 0 < hyph_row.angle < 90 or 180 < hyph_row.angle < 270:
        x0, y0 = hyph_row.x_mid - delta_x, hyph_row.y_mid - delta_y
        x1, y1 = hyph_row.x_mid + delta_x, hyph_row.y_mid + delta_y
    else:
        x0, y0 = hyph_row.x_mid - delta_x, hyph_row.y_mid + delta_y
        x1, y1 = hyph_row.x_mid + delta_x, hyph_row.y_mid - delta_y

    return x0, y0, x1, y1


def plot_to_pdf(snapshots, St_snapshots, S0, r, tstep, out_pdf):
    # dimensions for subplots on one page (n-rows and m-cols)
    n, m = 3, 4
    pdf_name = out_pdf
    with PdfPages(pdf_name) as pdf:
        # initialize layout for plots
        f, axarr = plt.subplots(n, m, sharex='none', sharey='none')
        arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
        subplots = [axarr[index] for index in arr_ij]

        splot_index = 0

        S_max_index = np.argmax(St_snapshots[max(St_snapshots)])

        for i in range(0, len(snapshots)):
            snapshot_to_plot = snapshots[sorted(snapshots)[i]]
            St_to_plot = St_snapshots[sorted(snapshots)[i]]

            endpoints_df = pd.DataFrame({'x0y0x1y1': snapshot_to_plot.apply(midpoints_to_endpoints, axis=1)})
            x0, y0, x1, y1 = endpoints_df['x0y0x1y1'].str[0].to_list(), endpoints_df['x0y0x1y1'].str[1].to_list(), \
                             endpoints_df['x0y0x1y1'].str[2].to_list(), endpoints_df['x0y0x1y1'].str[3].to_list()

            line_segments = [[list(zip(x0, y0))[i], list(zip(x1,y1))[i]] for i in
                             range(len(snapshot_to_plot))]
            # create continuous norm for mapping colors to hyphae according to time they occured
            hyphal_elements_norm = plt.Normalize(-10, len(snapshots)*1.25)
            red_cmap = clr.LinearSegmentedColormap.from_list('custom red', [(0, '#690909'), (0.5, '#AE1212'), (0.75, '#FD8D4C'),
                                                                (1, '#FCE0CB')],N=400)
            lc = mc.LineCollection(line_segments, linewidths=0.5, norm=hyphal_elements_norm, cmap=red_cmap) # OrRd_r
            lc.set_array(snapshot_to_plot['time'])
            subplots[splot_index].add_collection(lc)

            substrate_norm = plt.Normalize(0, S0)
            xlist = np.linspace(-20, 20, 100)
            ylist = np.linspace(-20, 20, 100)
            X, Y = np.meshgrid(xlist, ylist)
            dist = np.sqrt(X ** 2 + Y ** 2)
            dist[dist > 20] = 20
            Z = [[St_to_plot[int(r/20*x)] for x in y] for y in dist]
            subplots[splot_index].contourf(X, Y, Z, cmap=cm.get_cmap('YlGn'), norm=substrate_norm, alpha=0.7)

            subplots[splot_index].autoscale()
            subplots[splot_index].margins(0.1)
            subplots[splot_index].set_title(
                f"Time: {round(2 * i * tstep, 2)} h\n"
                f"Number of hyphal_elements: {len(snapshot_to_plot)}\n"
                f"Number of hyphal tips: {len(snapshot_to_plot[snapshot_to_plot.tip == True])}\n"
                f"Branching frequency: {round((len(snapshot_to_plot[snapshot_to_plot.tip == True]) / len(snapshot_to_plot)) / (tstep * i + float('1e-10')), 3)}\n"
                f"(branches per hyphal element per hour)",
                fontsize=5)
            subplots[splot_index].set_xlim(-20, 20)
            subplots[splot_index].set_ylim(-20, 20)
            subplots[splot_index].axis('off')
            #subplots[splot_index].set_xticks(list(range(-60,70,30)))
            #subplots[splot_index].set_yticks(list(range(-60,70,30)))

            plt.tight_layout()

            # increment index counter
            splot_index += 1
            # a pdf page is full once the count is equal to n*m
            if splot_index == n * m:
                pdf.savefig()
                plt.close(f)
                f, axarr = plt.subplots(n, m, sharex='none', sharey='none')
                arr_ij = [(x, y) for x, y in np.ndindex(axarr.shape)]
                subplots = [axarr[index] for index in arr_ij]
                splot_index = 0

        # save last page
        pdf.savefig()
        plt.close(f)
    return pdf_name


def plot_for_animation(snapshots, St_snapshots, S0, dirname, r, tstep, max_time):

    for i in range(len(snapshots)):
        snapshot_to_plot = snapshots[sorted(snapshots)[i]]
        St_to_plot = St_snapshots[sorted(snapshots)[i]]

        endpoints_df = pd.DataFrame({'x0y0x1y1': snapshot_to_plot.apply(midpoints_to_endpoints, axis=1)})
        x0, y0, x1, y1 = endpoints_df['x0y0x1y1'].str[0].to_list(), endpoints_df['x0y0x1y1'].str[1].to_list(), \
                         endpoints_df['x0y0x1y1'].str[2].to_list(), endpoints_df['x0y0x1y1'].str[3].to_list()

        fig, ax = plt.subplots(figsize=(8,8))

        line_segments = [[list(zip(x0, y0))[i], list(zip(x1, y1))[i]] for i in
                         range(len(snapshot_to_plot))]

        # create continuous norm for mapping colors to hyphae according to time they occured
        hyphal_elements_norm = plt.Normalize(0, max_time)
        purple_cmap = clr.LinearSegmentedColormap.from_list('custom purple', [(0, '#440154FF'),(0.75, '#404788FF'), (0.95, '#9999FF'), (1,'#CCCCFF')], N=200)
        lc = mc.LineCollection(line_segments, linewidths=1, norm=hyphal_elements_norm,
                               cmap=purple_cmap)  # cmap='OrRd_r'
        lc.set_array(snapshot_to_plot['time'])
        ax.add_collection(lc)

        substrate_norm = plt.Normalize(0, S0)  # S0
        green_cmap = clr.LinearSegmentedColormap.from_list('custom green', ['#FFFFFF', '#55C667FF'], N=200) #73D055FF
        xlist = np.linspace(-20, 20, 1000)
        ylist = np.linspace(-20, 20, 1000)
        X, Y = np.meshgrid(xlist, ylist)
        dist = np.sqrt(X ** 2 + Y ** 2)
        dist[dist>20] = 20
        Z = [[St_to_plot[int(r/20*x)] for x in y] for y in dist]

        ax.contourf(X, Y, Z, cmap=green_cmap, norm=substrate_norm,
                                       alpha=0.5)  # norm=substrate_norm    #YlGn cm.get_cmap('YlGn')

        ax.autoscale()
        #ax.margins(0.1)
        ax.set_title(
            f"Time: {round(i*tstep,2)} h\n"
            f"Number of hyphal_elements: {len(snapshot_to_plot)}\n"
            f"Number of hyphal tips: {len(snapshot_to_plot[snapshot_to_plot.tip == True])}\n"
            f"Branching frequency: {round((len(snapshot_to_plot[snapshot_to_plot.tip == True]) / len(snapshot_to_plot))/(tstep*i+float('1e-10')),3)}\n"
            f"(branches per hyphal element per hour)",
            fontsize=9)
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
        plt.axis('off')

        plt.tight_layout()
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        plt.savefig(f"{dirname}/image-{'0'*(4-len(str(i)))}{i}.png")
        plt.close(fig)

    return dirname


def make_gif(filenames, gifname):
    """
    Make a gif of files in filenames
    """
    with io.get_writer(gifname, mode='I', duration=0.05) as writer:
        for filename in filenames:
            image = io.imread(filename)
            writer.append_data(image)
    writer.close()
    return


def make_mp4(dirname):
    """
    Make an mp4 movie of png files in directory dirname
    """
    os.system(f"ffmpeg -framerate 20 -pattern_type glob -i '{dirname}/image*.png'"
              f" -c:v libx264 -r 30 -pix_fmt yuv420p {dirname}/movie.mp4")
    return


def main(args):

    # PARAMETERS USED IN THE SIMULATION

    # if given output from image analysis tool, use parameters from this
    if args.img_ana:
        img_ana_params = img_ana_parameters(args.img_ana)
        q = img_ana_params['branching_frequency'].iloc[0]
        curvature_gamma_params = (img_ana_params['curvature_gamma_a'].iloc[0], img_ana_params['curvature_gamma_scale'].iloc[0])
        branch_angle_beta_params = (img_ana_params['angle_beta_a'].iloc[0], img_ana_params['angle_beta_b'].iloc[0])

    # default distribution parameters are based on image analysis of ATCC 1015
    else:
        q = args.q
        curvature_gamma_params = (1.6114694580278444, 2.0238222211970487)
        branch_angle_beta_params = (4.28868341157166, 0.9839464337063648)


    # parameters for branching/extension
    p_lateral = 3/4*q                   # lateral branching frequency
    p_apical = 1/4*q                    # apical branching frequency
    lat_sub_min = args.lat_sub_min      # minimum substrate level for lateral branching
    ap_sub_min = args.ap_sub_min        # minimum substrate level for apical branching
    ext_sub_min = args.ext_sub_min       # minimum substrate level for tip extension
    mu_max = args.mu_max                # maximal growth rate (exponential phase)

    # parameters for substrate gradient
    fieldtype = args.fieldtype          # type of substrate field (uniform or gradient)
    source = args.source                # source / no source at edge of field
    S0 = args.S0                        # initial substrate concentration
    r = 400                             # "radius" of substrate field (2000/100 = 20)
    h = 0.01                            # radial step size
    D = args.D                          # diffusion constant
    S_tip = args.S_tip                  # substrate consumption by tips
    S_nontip = args.S_nontip            # substrate consumption by non-tips
    tstep = args.tstep                  # time of 1 simulation round in hours
    N = args.hours/tstep                # max simulation rounds
    branch_substrate_dependency = 1.3   # branching dependency on substrate concentration
    St = [S0 for i in range(r + 1)]     # initial uniform substrate field


    # lists for saving results during the simulation
    snapshots = dict()
    St_snapshots = dict()

    # initialize the first two hyphal_elements growing from a single spore
    n_hyphal_elements = 8
    hyphal_elements_dict = {'x_mid':[-0.5,0.5,0,0, np.sqrt(2)/4, -np.sqrt(2)/4, -np.sqrt(2)/4, np.sqrt(2)/4], 'y_mid':[0,0,-0.5,0.5, np.sqrt(2)/4, np.sqrt(2)/4, -np.sqrt(2)/4, -np.sqrt(2)/4], 'angle':[180,0,270,90,45,135,225,315 ], 'tip':[True, True, True, True, True, True, True, True], 'time':[0,0,0,0,0,0,0,0]}
    hyphal_elements = pd.DataFrame(hyphal_elements_dict)

    # append initial dataframe to snapshot lists
    snapshots[0] = (copy.deepcopy(hyphal_elements))
    St_snapshots[0] = St

    # actual simulation
    i, m = 1, 0
    while i < N:
        hyphal_elements_tip = hyphal_elements[hyphal_elements.tip == True]
        hyphal_elements_nontip = hyphal_elements[hyphal_elements.tip == False]
        St = substrate_per_distance(n_tip=len(hyphal_elements_tip), n_nontip=len(hyphal_elements_nontip), St=St, D=D, S_tip=S_tip, S_nontip=S_nontip, r=r, fieldtype = fieldtype, source = source)

        for j in range(len(hyphal_elements_nontip)):
            # lateral branching for non-tip hyphal_elements
            branching(hyphal_elements_nontip.iloc[j], St, lat_sub_min, (branch_substrate_dependency*np.mean(St)/S0)*p_lateral, hyphal_elements, branch_angle_beta_params, r=r, time=i)

        for k in range(len(hyphal_elements_tip)):
            # apical branching for hyphal tips
            branching(hyphal_elements_tip.iloc[k], St, ap_sub_min, p_apical, hyphal_elements, branch_angle_beta_params, r=r, time=i)
            # extensions for hyphal tips, using the Monod equation
            tip_extension_monod(mu_max, hyph_row=hyphal_elements_tip.iloc[k], hyphal_elements=hyphal_elements, St=St, Ks=200,
                                ext_sub_min=ext_sub_min, curvature_gamma_params=curvature_gamma_params, r=r, time=i)

        i += 1
        if (i) % 2 == 0:
            print(f"Iteration = {i}\tNumber of hyphal_elements = {len(hyphal_elements)}")
            snapshots[i] = copy.deepcopy(hyphal_elements)
            St_snapshots[i] = copy.deepcopy(St)

    print('# Simulation done! Plotting...')


    # make gif animation
    if args.ani:
        dir_name = plot_for_animation(snapshots, St_snapshots, S0, dirname=args.ani, r=r, tstep=tstep, max_time=i)
        filenames = sorted(dir_name+'/'+fn for fn in os.listdir(dir_name) if fn.startswith('image'))
        make_gif(filenames, f'{dir_name}/{args.ani}.gif')
        make_mp4(dir_name)

    # make pdf of subplots
    if args.pdf:
        pdf_name = plot_to_pdf(snapshots, St_snapshots, S0, r=r, tstep=tstep, out_pdf=args.pdf)

    return


if __name__ == '__main__':
    print("# Running Mycemulator...")
    start_time = datetime.now()
    args = get_args()
    print("# arguments:", args)
    pdf_name = main(args)
    end_time = datetime.now()
    print(f'# Done!\nDuration: {end_time-start_time}')
#!/usr/bin/python
"""
Plot network evaluations in rank-precision style. 
Example usage:
python plot_evaluations.py --network_evals <network_eval_list> --network_labels <label_list> --random_eval_dir <random_eval_dir> --figure_file_suffix <figure_suffix> --step 1600 --num_regulators 320
"""

import sys
import numpy as np
import argparse
import glob
import os.path
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


EVAL_METHODS = ['binned', 'cumulative']
## color theme
COLOR_THEME = {"black":(0,0,0), 
        "blue":(31, 119, 180), "blue_L":(174, 199, 232), 
        "orange":(255, 127, 14), "orange_L":(255, 187, 120),
        "green":(44, 160, 44), "green_L":(152, 223, 138), 
        "red":(214, 39, 40), "red_L":(255, 152, 150),
        "magenta":(148, 103, 189), "magenta_L":(197, 176, 213),
        "brown":(140, 86, 75), "brown_L":(196, 156, 148),
        "pink":(227, 119, 194), "pink_L":(247, 182, 210), 
        "grey":(127, 127, 127), "grey_L":(199, 199, 199),
        "yellow":(255, 215, 0), "yellow_L":(219, 219, 141), 
        "cyan":(23, 190, 207), "cyan_L":(158, 218, 229)}
ORDERED_COLORS = ["blue", "orange", "green", "red", "magenta", 
                  "brown", "pink", "grey", "yellow", "cyan"]
for c in COLOR_THEME.keys():
    (r, g, b) = COLOR_THEME[c]
    COLOR_THEME[c] = (r/255., g/255., b/255.)


def parse_args(argv):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--network_evals', 
                        help="List of network evaluation files, delimited by ','")
    parser.add_argument('--network_labels', 
                        help="List of network labels, delimited by ','")
    parser.add_argument('--random_eval_dir', default=None,
                        help="Directory of evaluations of randominzed networks")
    parser.add_argument('--figure_file_suffix', 
                        help="Suffix of figure file. Two files will be created: ChIP support and PWM support.")
    parser.add_argument('--step', default=1600, type=float)
    parser.add_argument('--num_regulators', default=320, type=int)
    parser.add_argument('--eval_method', default='cumulative', help='%s' % EVAL_METHODS)
    parsed = parser.parse_args(argv[1:])
    return parsed


def main(argv):
    parsed = parse_args(argv)
    
    ## Create evaluation of true networks
    network_evals = [x.strip() for x in parsed.network_evals.split(",")]
    network_labels = [x.strip() for x in parsed.network_labels.split(",")]
    if len(network_evals) > len(ORDERED_COLORS):
        sys.exit("Too many networks for color handling. Abort.")
    network_colors = ORDERED_COLORS[:len(network_evals)]

    ## Create chance evaluation of randomized networks
    if parsed.random_eval_dir is not None:
        rand_eval_chip = []
        rand_eval_pwm = []
        for file_rand in sorted(glob.glob(parsed.random_eval_dir + "/*")):
            chip_out, pwm_out = parse_binding_overlap(file_rand, parsed.eval_method)
            rand_eval_chip.append(chip_out)
            rand_eval_pwm.append(pwm_out)
        chance_eval_chip = [numpy.percentile(rand_eval_chip, 2.5, axis=0),
            numpy.percentile(numpy.array(rand_eval_chip), 97.5, axis=0)]
        chance_eval_pwm = [numpy.percentile(rand_eval_pwm, 2.5, axis=0), 
            numpy.percentile(numpy.array(rand_eval_pwm), 97.5, axis=0)]
    else:
        chance_eval_chip, chance_eval_pwm = None, None

    ## Make rank precision plots
    plot_analysis(network_evals, 
                network_colors, 
                network_labels, 
                parsed.figure_file_suffix, 
                parsed.num_regulators, 
                parsed.step, 
                parsed.eval_method, 
                chance_eval_chip, 
                chance_eval_pwm)


def plot_analysis(network_evals, colors, labels, figure_file_suffix, num_regulators, step, eval_method, chance_eval_chip, chance_eval_pwm):
    # compute chip and pwm supports
    eval_chip = [None] * len(network_evals)
    eval_pwm = [None] * len(network_evals)
    exp_chance_chip, exp_chance_pwm = parse_chance_binding_overlap(network_evals[0])
    for i in range(len(network_evals)):
        [eval_chip[i], eval_pwm[i]] = parse_binding_overlap(network_evals[i], eval_method)
    
    x_ticks = [format(float(i)*step/num_regulators, '.0f') for i in range(1,len(eval_chip[0])+1)]  

    font_size = 14
    matplotlib.rcParams.update({'font.size': 16})  

    ## ChIP support
    fig, ax = plt.subplots(figsize=(7,6), dpi=150)
    for i in range(len(eval_chip)):
        ax.plot(eval_chip[i], color=colors[i], label=labels[i], linewidth=3)
    if chance_eval_chip is not None:
        ax.fill_between(numpy.arange(20), chance_eval_chip[0], chance_eval_chip[1],
                    facecolor="black", alpha=.25, label="Random: 95% CI")
    ax.plot(exp_chance_chip, color="black", linestyle=":", linewidth=3,
                    label="Random: expectation",)
    
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by ChIP (%)')
    ax.yaxis.grid(True, linestyle='-', alpha=.2, linewidth=.5)
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, 50)
    plt.yticks(numpy.arange(0,50.5,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(fontsize=font_size)
    plt.savefig(figure_file_suffix + '.ChIP_support.pdf', fmt='pdf')

    ## PWM support
    fig, ax = plt.subplots(figsize=(7,6), dpi=150)
    for i in range(len(eval_pwm)):
        ax.plot(eval_pwm[i], color=colors[i], label=labels[i], linewidth=3)
    if chance_eval_pwm is not None:
        ax.fill_between(numpy.arange(20), chance_eval_pwm[0], chance_eval_pwm[1], 
                    facecolor="black", alpha=.25, label="Random: 95% CI")
    ax.plot(exp_chance_pwm, color="black", linestyle=":", linewidth=3, 
                    label="Random: expectation")

    ax.scatter(18, 10, s=75, c=COLOR_THEME["yellow"])
    ax.annotate('ChIP network', xy=(17.95, 9.75), xycoords='data', xytext=(-50, -30), textcoords='offset points', arrowprops=dict(arrowstyle="->"))
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Average number of predicted targets per TF in the genome')
    plt.ylabel('Interactions supported by PWM (%)')
    ax.yaxis.grid(True, linestyle='-', alpha=.2, linewidth=.5)
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, 40)
    plt.yticks(numpy.arange(0,40.5,5))
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(fontsize=font_size)
    plt.savefig(figure_file_suffix + '.PWM_support.pdf', fmt='pdf')
    

def plot_bar_analysis(network_evals, colors, labels, figure_file_suffix, num_regulators, step, eval_method, chance_eval_chip, chance_eval_pwm):
    # compute chip and pwm supports
    eval_chip = [None] * (len(network_evals)+1)
    eval_pwm = [None] * (len(network_evals)+1)
    [eval_chip[0], eval_pwm[0]] = parse_chance_binding_overlap(network_evals[0])
    for i in range(len(network_evals)):
        [eval_chip[i+1], eval_pwm[i+1]] = parse_binding_overlap(network_evals[i], eval_method)

    # make bar plots
    indx = range(2)
    width = .35
    colors[0] = COLOR_THEME["grey"]
    fig, ax = plt.subplots(figsize=(5,4.5), dpi=150)
    for i in range(1,len(eval_chip)):
        ax.bar([ind*2+width*(i+2) for ind in indx], (eval_chip[i][1], eval_pwm[i][1]), width, color=colors[i], label=labels[i])
    ax.set_ylabel("Interactions supported by benchmarks (%)")
    ax.set_xticks([ind*2+width*(i+.5) for ind in indx])
    ax.set_xticklabels(("ChIP", "PWM"))
    ax.set_ylim([0, 40])
    ax.legend(prop={'size':12})
    plt.savefig(figure_file_suffix + '.bar_top_3200.pdf', fmt='pdf')


def parse_binary_gold_standard(filenames, method):
    eval_chip = numpy.zeros([len(filenames)/2+1, 10])
    eval_pwm = numpy.zeros([len(filenames)/2+1, 10])

    for i in range(len(filenames)/2):
        chip_support = numpy.loadtxt(filenames[i*2])
        pwm_support = numpy.loadtxt(filenames[i*2+1]) 
        if i == 0:
            eval_chip[0,:] = chip_support[0,:]
            eval_pwm[0,:] = pwm_support[0,:]
        eval_chip[i+1,:] = chip_support[1,:]
        eval_pwm[i+1,:] = pwm_support[1,:]
        
    if method == 'cumulative':
        temp_eval_chip = numpy.zeros([len(filenames)/2+1, 10])
        temp_eval_pwm = numpy.zeros([len(filenames)/2+1, 10])
        for j in range(10):
            temp_eval_chip[:,j] = numpy.sum(eval_chip[:,0:(j+1)], axis=1)/(j+1)
            temp_eval_pwm[:,j] = numpy.sum(eval_pwm[:,0:(j+1)], axis=1)/(j+1)
        eval_chip = temp_eval_chip
        eval_pwm = temp_eval_pwm

    return [eval_chip*100, eval_pwm*100]


def parse_binding_overlap(filename, method):
    lines = open(filename, "r").readlines()
    chip = [0] * (len(lines))
    pwm = [0] * (len(lines))

    if method == "cumulative":
        for i in range(len(lines)):
            line = lines[i].split()
            chip[i] = float(line[5])/float(line[2])
            pwm[i] = float(line[4])/float(line[2])
    
    elif method == "binned":
        for i in range(len(lines)):
            line = lines[i].split()
            if i == 0:
                chip[i] = float(line[5])/float(line[2])
                pwm[i] = float(line[4])/float(line[2])
            else:
                chip[i] = (float(line[5]) - float(prevline[5]))/(float(line[2]) - float(prevline[2]))
                pwm[i] = (float(line[4]) - float(prevline[4]))/(float(line[2]) - float(prevline[2]))
            prevline = line  
    return [numpy.array(chip)*100, numpy.array(pwm)*100]


def parse_chance_binding_overlap(fn):
    evalPoints = numpy.loadtxt(fn).shape[0]
    line = open(fn, 'r').readline()
    line = line.split()
    chip = [float(line[1])/float(line[7]) for _ in range(evalPoints)]
    pwm = [float(line[0])/float(line[7]) for _ in range(evalPoints)]
    return [numpy.array(chip)*100, numpy.array(pwm)*100]


if __name__ == "__main__":
    main(sys.argv)

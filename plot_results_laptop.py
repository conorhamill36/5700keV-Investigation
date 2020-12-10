#python script to plot results from output file of state_finder.C, convert angle into CoM frame, and plot alongside results from fresco calculations
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
matplotlib.rcParams.update({'errorbar.capsize': 2}) #Copied from  https://stackoverflow.com/questions/18434492/matplotlib-errorbar-caps-missing to fix error bar caps not appearing


def plot_energy_diff(energy_diff_angle_array, energy_diff_com_angle_array, energy_diff_array, energy_diff_uncert_array):
    print("Plotting energy diff ")
    print("Trimming arrays")
    energy_diff_angle_array = energy_diff_angle_array[3:-2]
    energy_diff_com_angle_array = energy_diff_com_angle_array[3:-2]
    energy_diff_array = energy_diff_array[3:-2]
    energy_diff_uncert_array = energy_diff_uncert_array[3:-2]


    plt.errorbar(energy_diff_com_angle_array, energy_diff_array, energy_diff_uncert_array, \
    color='firebrick', barsabove=True, marker='.', ms=7, linestyle='None', elinewidth=1, capthick=1,\
    markeredgewidth=2)
    plt.ylabel(r'$\mathrm{\Delta}$E [keV]', fontsize=15,rotation=90, labelpad=-2)
    plt.xlabel(r'$\mathrm{\theta_{CoM}[^{\circ}]}$', fontsize=15, labelpad=-2)
    plt.xlim(0,60)

    plt.gca().set_xticks(np.arange(0, 60, 10))

    # plt.xticks([])
    plt.xlabel(r'$\mathrm{\theta_{CoM}[^{\circ}]}$', fontsize=15, labelpad=-2)
    # fig.subplots_adjust(wspace=None, hspace=None)
    # plt.tight_layout()
    # plt.show()
    return

def draw_arrows():


    #Above plot
    # l_1_peak, l_2_peak, l_3_peak = 16, 30, 45
    # arrow_base_y = 27.2
    # arrow_height = -0.7
    # head_width = -arrow_height*1.5
    # head_length = head_width * 0.5

    #Below plot
    l_1_peak, l_2_peak, l_3_peak = 16, 30, 45
    arrow_base_y = 19.2
    arrow_height = 0.4
    head_width = arrow_height*1.5
    head_length = head_width * 0.5

    #l=1 arrow
    plt.arrow(l_1_peak, arrow_base_y, 0, arrow_height, clip_on = False, shape='full',\
    head_width = head_width, head_length = head_length)
    #l=1 label
    plt.text(l_1_peak, arrow_base_y + 0.2, r'$\ell=1$ peak', ha = 'center', wrap=True)

    #l=2 arrow
    plt.arrow(l_2_peak, arrow_base_y, 0, arrow_height, clip_on = False, shape='full',\
    head_width = head_width, head_length = head_length)
    #l=2 label
    plt.text(l_2_peak, arrow_base_y + 0.2, r'$\ell=2$ peak', ha = 'center')



    #l=3 arrow
    plt.arrow(l_3_peak, arrow_base_y, 0, arrow_height, clip_on = False, shape='full',\
    head_width = head_width, head_length = head_length)
    #l=3 label
    plt.text(l_3_peak, arrow_base_y + 0.2, r'$\ell=3$ peak', ha = 'center')


    return




def main():

    #Opening input file
    plot_results_input_file = open("plot_results_input.txt","r")
    plot_results_input = plot_results_input_file.readlines()
    i=0


    plot_results_input_file = open("exp_fresco_diff.txt","w+") #nasty hack to empty file every time
    plot_results_input_file.close()
    # fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True)



    for line in plot_results_input:

        # ax = axs[i]
        exp_fresco_diff_output = open("exp_fresco_diff.txt","a")
        if( (line.startswith('Ex') == 0)  and (line.startswith('#') == 0)):
            token=line.split("\t")
            print(token)
            E_x = token[0]
            J_pi = token[1]
            L_zero = float(token[2])
            L_one = float(token[3])
            L_two = float(token[4])
            L_three = float(token[5])
            S_0 = float(token[6])
            S_1 = float(token[7])
            S_2 = float(token[8])
            S_3 = float(token[9])
            state_compound_factor = float(token[10])
            ax = plot_single_state(E_x, J_pi, L_zero, L_one, L_two, L_three, S_0, S_1, S_2, S_3, i, state_compound_factor)
            i=i+1
    plot_results_input_file.close()

    lab_cross_section = 1.0
    E_ex = 5.0
    CoM_angle = 10

    CoM_cross_section = cross_section_lab_to_CoM_frame(lab_cross_section, E_ex, CoM_angle)

    print(CoM_cross_section)
    exp_fresco_diff_output.close()

def cross_section_lab_to_CoM_frame(lab_cross_section, E_ex, CoM_angle):
    ma = 2.01410177811
    mA = 24.98583696
    mb = 1.00782503224
    mB = 25.98259297
    E_beam = 8.0
    Q_gs = 8.868514

    Q = Q_gs - float(E_ex)
    # print("Q is {}".format(Q))

    gamma = math.sqrt((ma*mb*E_beam) / (mB * (mb + mB)*Q + mB*(mB+mb-ma)*E_beam))
    transform_factor = pow(1 + gamma*gamma + 2*gamma*math.cos(math.radians(CoM_angle)), 1.5) / (1 + gamma*(math.cos(math.radians(CoM_angle))))
    CoM_cross_section = lab_cross_section * transform_factor
    # print("gamma is {}, transform factor is {}, gicing a com x s of {}".format(gamma,transform_factor,CoM_cross_section))
    return CoM_cross_section

def plot_single_state(E_x,J_pi,L_zero,L_one,L_two,L_three,S_0,S_1,S_2,S_3,counter, state_compound_factor):

    #Setting arrays
    lab_angles_list=[]
    com_angles_list=[]
    target_4_lab_angles_list=[]
    target_4_cs_list=[]
    target_4_com_angles_list=[]
    target_4_cs_uncert_list=[]

    target_5_lab_angles_list=[]
    target_5_com_angles_list=[]
    target_5_cs_list=[]
    target_5_cs_uncert_list=[]
    targets_com_angles_list = []
    targets_cs_list = []
    targets_cs_uncertainty_list = []
    compound_angles_list = []
    compound_cs_list = []
    compound_fresco_cs_list = []
    temp_compound_fresco_cs_list = [0 for i in range (181)]
    OMP_angle_list = []
    OMP_stddev_list =[]
    OMP_angle_L_0_list = []
    OMP_stddev_L_0_list = []
    OMP_angle_L_2_list = []
    OMP_stddev_L_2_list = []
    yasue_com_angle_list = []
    yasue_cs_list = []

    n=4
    m=180
    fresco_angle_list=[[] * m for i in range(n)]
    fresco_cross_section_list=[[] * m for i in range(n)]
    fresco_angle_sum_list=[0 for i in range (181)]
    fresco_cross_section_sum_list=[0 for i in range (181)]
    upper_spec_list=[0 for i in range (181)]
    lower_spec_list=[0 for i in range (181)]
    exp_fresco_diff_list = [7]
    exp_fresco_diff_list_angle = [7]

    temp_fresco_angle_list=[[] * m for i in range(n)]
    temp_fresco_cross_section_list=[[] * m for i in range(n)]
    temp_fresco_angle_sum_list=[0 for i in range (181)]

    #variables
    target_4_flag=0
    target_5_flag=0

    #booleans
    save_fig = 0 #saves figure to pdf and won't show on screen
    combined_fit = 0 #saves in different name - need to check
    convert_to_CoM_frame = 1 #converts cross sections to CoM frame
    compound_boolean = 1 #includes compound nuclear component of reaction
    add_compound_boolean = 1 #adds compound nuclear component of reaction to total model cross-section
    fit_boolean = 1 #fits spectroscopic factor to experimental data
    diff_boolean = 0 # calculates difference between experiment and fresco
    target_thickness_uncertainty_boolean = 0 #Adds 10% uncertainty from target thickness to plot
    OMP_uncertainty_boolean = 0
    chi_square_boolean = 1 #Calculates spec factor 1 sigma outside chi-square minimum
    yasue_boolean = 1 #plots Yasue's 1+ 5.691 MeV spec factor of S=0.2 to plot.
    paper_boolean = 0 #turns on and off aesthetic things to get plot ready for paper
    fitting_uncertainty_boolean = 1 #turns on and off adding uncertainty from choosing fitting parameters
    compound_nuclear_factor = 0.10 #factor compound nuclear component from TALYS output is multiplied by before being used
    positive_parity_boolean = 1
    negative_parity_boolean = 0
    both_parities_boolean = 1 #boolean to investigate negative contributions for 5.7 MeV state
    #Manual inputs
    #allowed angular momentum contributions
    # L_zero = 0
    # L_one = 1
    # L_two = 0
    # L_three = 1
    #
    #
    #Spectroscopic Factors Used (Values of S, not (2J+1)S)
    # S_0 = 0.0
    # S_1 = 0.095
    # S_2 = 0.0
    # S_3 = 0.04
    #
    #
    # E_x = "6.876" #Excitations energy in MeV
    # J_pi = "3-" #Angular momentum and parity of state
    # delta_L = "1"

    L_list=[L_zero, L_one, L_two, L_three]
    spec_list=[S_0, S_1, S_2, S_3]

    J = float(J_pi[:1])
    print('J is {}'.format(J))



    #Transform to CoM angles
    #Using table with transformation for Ex=6.0 MeV
    angles_file = open("com_angle.txt","r")
    angles_input = angles_file.readlines()
    exp_fresco_diff_output = open("exp_fresco_diff.txt","a")
    for line in angles_input:
        # print(line)
        token=line.split("\t")
        # print(token[0])
        token[0]=float(token[0])
        # print(token[1])
        token[1]=float(token[1])
        lab_angles_list.append(token[0])
        com_angles_list.append(token[1])
    # print(lab_angles_list)
    # print(com_angles_list)

    #result_file_name = "results/" + E_x + "_" + J_pi + "_results.txt"
    result_file_name = "5.716_4+_results.txt"
    print("Results file name is {}".format(result_file_name))
    if(combined_fit):
        result_file_name = "results/" + E_x + "_" + J_pi + "_combined_results.txt"

    # result_file = open("results/6.876_3-_results.txt")
    result_file = open(result_file_name)
    #Taking experimental cross sections from file
    results_input = result_file.readlines()
    for line in results_input:
        # print(line)
        if line.startswith("Target#4"):
            print("Target 4 found, flag set")
            target_4_flag=1
            target_5_flag=0
            continue
        if line.startswith("Target#5"):
            print("Target 5 found, flag set")
            target_4_flag=0
            target_5_flag=1
            continue

        if line.startswith("0 "):
            continue

        token=line.split(" ")
        token[0]=float(token[0])
        token[1]=float(token[1])
        token[2]=float(token[2])

        if(target_4_flag==1):
            target_4_lab_angles_list.append(token[0])
            target_4_cs_list.append(token[1])
            target_4_cs_uncert_list.append(token[2])

        if(target_5_flag==1):
            target_5_cs_list.append(token[1])
            target_5_lab_angles_list.append(token[0])
            target_5_cs_uncert_list.append(token[2])



    #Removing first three angles we don't want to plot
    for i in range (0,4):
        # print(target_4_lab_angles_list[i])
        dummy_variable=target_4_lab_angles_list[0]
        if (dummy_variable==5.0 or dummy_variable==7.0 or dummy_variable==10.0):
            target_4_cs_list.pop(0)
            target_4_lab_angles_list.pop(0)
            target_4_cs_uncert_list.pop(0)
            # target_4_cs_list[i]=-10.0
    #Removing points for merged angles
    if(E_x == '6.978'): #Removes values at 13, 16, 19 degrees
        for m in range(0,3):
            target_4_cs_list.pop(0)
            target_4_lab_angles_list.pop(0)
            target_4_cs_uncert_list.pop(0)

    if(E_x == '7.061'): #Removes values at 30, 33, 36, 39 degrees
        target_4_cs_list.pop(6)
        target_4_lab_angles_list.pop(6)
        target_4_cs_uncert_list.pop(6)
        for m in range(0,4):
            target_5_cs_list.pop(0)
            target_5_lab_angles_list.pop(0)
            target_5_cs_uncert_list.pop(0)
        target_5_cs_list.pop(2)
        target_5_lab_angles_list.pop(2)
        target_5_cs_uncert_list.pop(2)

    if(E_x == '7.349'): #Removes values at 36, 39, 45,
        for m in range(0,4):
            target_5_cs_list.pop(2)
            target_5_lab_angles_list.pop(2)
            target_5_cs_uncert_list.pop(2)
    if(E_x == '7.100'):
        for m in range(0,6):
            target_5_cs_list.pop(2)
            target_5_lab_angles_list.pop(2)
            target_5_cs_uncert_list.pop(2)

    if(E_x == '7.261'):
        for m in range(0,2):
            target_5_cs_list.pop(4)
            target_5_lab_angles_list.pop(4)
            target_5_cs_uncert_list.pop(4)
    if(E_x == '5.292'):
            target_4_cs_list.pop(4)
            target_4_lab_angles_list.pop(4)
            target_4_cs_uncert_list.pop(4)
    if(E_x == '5.716'):
        target_4_cs_list.pop(3)
        target_4_lab_angles_list.pop(3)
        target_4_cs_uncert_list.pop(3)



    print(target_4_lab_angles_list)
    print(target_4_cs_list)

    print(target_5_lab_angles_list)
    print(target_5_cs_list)

    angles_file.close()
    result_file.close()

    #Converting lab angles to CoM angles
    #Target 4 conversion
    print("target 4 conversion")
    # print(len(target_4_lab_angles_list))
    for i in range (0,len(target_4_lab_angles_list)):
        # print(target_4_lab_angles_list[i])
        if target_4_lab_angles_list[i] in lab_angles_list:
            # print("match found!")
            com_index = lab_angles_list.index(target_4_lab_angles_list[i])
            # print(com_index)
            target_4_com_angles_list.append(com_angles_list[com_index])

    print("target 5 conversion")
    # print(len(target_5_lab_angles_list))
    for i in range (0,len(target_5_lab_angles_list)):
        # print(target_5_lab_angles_list[i])
        if target_5_lab_angles_list[i] in lab_angles_list:
            # print("match found!")
            com_index = lab_angles_list.index(target_5_lab_angles_list[i])
            # print(com_index)
            target_5_com_angles_list.append(com_angles_list[com_index])

    # print('\n')
    # # print(target_4_lab_angles_list)
    # # print(target_4_com_angles_list)
    # print('\n')
    #
    #
    # print('\n')
    # # print(target_5_lab_angles_list)
    # # print(target_5_com_angles_list)
    print('\n')

    if(convert_to_CoM_frame):
        #Converting cross-sections from lab-frame to CoM frame
        print("Converting to CoM frame")
        for i in range(0, len(target_4_com_angles_list)):
            target_4_cs_list[i] = cross_section_lab_to_CoM_frame(target_4_cs_list[i], E_x, target_4_com_angles_list[i])
            print("{}\t{}".format(target_4_com_angles_list[i],target_4_cs_list[i]))

        for i in range(0, len(target_5_com_angles_list)):
            target_5_cs_list[i] = cross_section_lab_to_CoM_frame(target_5_cs_list[i], E_x, target_5_com_angles_list[i])
            print("{}\t{}".format(target_5_com_angles_list[i],target_5_cs_list[i]))

    if(target_thickness_uncertainty_boolean):
        target_thickness_uncert = 0.10 #10% uncertainty from target thickness
        print("adding target thickness of {}".format(target_thickness_uncert))
        for i in range (0,len(target_4_cs_uncert_list)):
            stat_uncert = target_4_cs_uncert_list[i]
            total_uncert = stat_uncert*stat_uncert + target_thickness_uncert*target_4_cs_list[i]*target_thickness_uncert*target_4_cs_list[i]
            total_uncert = math.sqrt(total_uncert)
            # print("stat uncert is {}, total uncert is now {}".format(stat_uncert, total_uncert))
            target_4_cs_uncert_list[i] = total_uncert


        for i in range (0,len(target_5_cs_uncert_list)):
            stat_uncert = target_5_cs_uncert_list[i]
            total_uncert = stat_uncert*stat_uncert + target_thickness_uncert*target_5_cs_list[i]*target_thickness_uncert*target_5_cs_list[i]
            total_uncert = math.sqrt(total_uncert)
            # print("stat uncert is {}, total uncert is now {}".format(stat_uncert, total_uncert))
            target_5_cs_uncert_list[i] = total_uncert




    #Finding fresco outputs and making lists for individual and sum contributions
    for i in range (0,4):
        if L_list[i] == 1:
            print("contribution from {} needed".format(i))
            delta_L=i
            delta_L=str(delta_L)
            fresco_file_name = "mg26dp_dwba_5.716_" + J_pi + "_" + delta_L + ".dat"
            #fresco_file_name="/home/s1771161/Mg25_d_p_Mg26/src/fresco_outputs_bcfec/mg26dp_dwba_" + E_x + "_" + J_pi + "_" + delta_L +".dat"
            fresco_file = open(fresco_file_name)
            fresco_data = fresco_file.readlines()
            j=0
            for line in fresco_data:
                token=line.split("\t")
                token[0]=float(token[0])
                token[1]=float(token[1])


                #Multiplying by spec factors

                if (spec_list[i] > 0.0001):
                    token[1]=token[1]*spec_list[i]

                fresco_angle_list[i].append(token[0])
                fresco_cross_section_list[i].append(token[1]) #i is L contribution index
                fresco_cross_section_sum_list[j]=fresco_cross_section_sum_list[j]+token[1] #j is angle index
                # fresco_angle_sum_list.append(token[0])
                j=j+1


            fresco_file.close()
            # print(fresco_angle_list[i])
            # print(fresco_cross_section_list[i])



    for i in range (0,181):
        fresco_angle_sum_list[i]=i


    #Finding compound nuclear cross sections and adding to list
    if(compound_boolean):
        #compound_file_name = "/scratch/TALYS/talys/mg25dpmg26/talys_outputs/" + E_x + "_talys.dat"
        compound_file_name = "5.716_talys.dat"
        print("compound cs's read in from {}".format(compound_file_name))
        compound_file = open(compound_file_name,"r")
        compound_input = compound_file.readlines()
        for line in compound_input:
            # print(line)
            token=line.split("\t")
            token[-1] = token[-1].strip()
            # print(token)
            compound_angles_list.append(float(token[0]))
            compound_cs_list.append(compound_nuclear_factor*state_compound_factor*float(token[1]))

        compound_file.close()

    if(add_compound_boolean): #adding compound and fresco sum together
        for m in range(0,176):
            # print(len(compound_fresco_cs_list))
            # compound_fresco_cs_list.append(compound_cs_list[m] + fresco_cross_section_list[2][m])
            compound_fresco_cs_list.append(compound_cs_list[m] + fresco_cross_section_sum_list[m])

    #Fitting experimental data to best spectroscopic factor
    # for k in range(0,50):
    #
    #     #Fitting L=0 contribution to 13 degrees point
    #     if(fit_boolean):
    #         if(spec_list[0]>0.0001):
    #             print("Fitting L=0 contribution")
    #             #Set new limits for spec factors
    #
    #             upper_spec = 1.4*spec_list[0]
    #             lower_spec = 0.6*spec_list[0]
    #
    #             for j in range(0,5): #loop that fits every time
    #                 print("\n\nlower spec is {}, upper spec is {}".format(lower_spec, upper_spec))
    #     #
    #                 # upper_spec_list[181] = 3.0*fresco_cross_section_list[2]
    #                 # lower_spec_list[181] = lower_spec*fresco_cross_section_list[2]
    #                 for i in range (0,175):
    #                     # print(fresco_cross_section_list[2][i])
    #                     # print(upper_spec_list[i])
    #                     if(add_compound_boolean):
    #                         upper_spec_list[i] = (upper_spec/spec_list[0])*compound_fresco_cs_list[i]
    #                         lower_spec_list[i] = (lower_spec/spec_list[0])*compound_fresco_cs_list[i]
    #                     else:
    #                         upper_spec_list[i] = (upper_spec/spec_list[0])*fresco_cross_section_list[0][i]
    #                         lower_spec_list[i] = (lower_spec/spec_list[0])*fresco_cross_section_list[0][i]
    #
    #     #             #Calculate difference between fresco and exp data points for both spec factor points
    #                 sum_diff = 0
    #                 #contribution from target 4
    #                 rounded_angle = round(target_4_com_angles_list[0])
    #                 # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
    #                 rounded_angle = int(rounded_angle)
    #                 print(rounded_angle)
    #                 # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
    #                 # if(add_compound_boolean):
    #                 # diff = target_4_cs_list[i] -
    #                 diff = target_4_cs_list[0] - upper_spec_list[rounded_angle]
    #                 diff = abs(diff)
    #                 sum_diff = sum_diff + abs(diff)
    #                     # print("diff is {}, giving sum diff of {}".format(diff, sum_diff))
    #
    #                 upper_sum_diff = sum_diff
    #
    #                 sum_diff = 0
    #                 #contribution from target 4
    #                 rounded_angle = round(target_4_com_angles_list[0])
    #                 # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
    #                 rounded_angle = int(rounded_angle)
    #                 # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
    #                 # if(add_compound_boolean):
    #                 # diff = target_4_cs_list[i] -
    #                 diff = target_4_cs_list[0] - lower_spec_list[rounded_angle]
    #                 diff = abs(diff)
    #                 sum_diff = sum_diff + abs(diff)
    #
    #                 lower_sum_diff = sum_diff
    #
    #                 print("lower sum diff is {}, upper sum diff is {}".format(lower_sum_diff,upper_sum_diff))
    #
    #                 #Decide which is closer and reset spec factor limits
    #                 if(lower_sum_diff < upper_sum_diff): #lower gives a better estimate
    #                     print("lower gives a better estimate, with a spec factor of {}".format(lower_spec))
    #                     lower_spec = lower_spec
    #                     # upper_spec = (upper_spec + lower_spec)/2.0
    #                     upper_spec = lower_spec + (upper_spec-lower_spec)*0.7
    #                     current_spec = lower_spec
    #
    #                 if(lower_sum_diff > upper_sum_diff): #upper gives a better estimate
    #                     print("upper gives a better estimate, with a sepc factor of {}".format(upper_spec))
    #                     upper_spec = upper_spec
    #                     # lower_spec = (lower_spec + upper_spec)/2.0
    #                     lower_spec = upper_spec - (upper_spec-lower_spec)*0.7
    #                     current_spec = upper_spec
    #
    #     #         #Outside fitting loop
    #             print("best spectroscopic factor is {}, input spectroscopic factor was {}".format(current_spec, spec_list[0]))
    #             S_0 = current_spec
    #     #         #changing spectroscopic factor that fresco output is scaled by
    #     #
    #             for i in range (0,174): #changing model output to new spectroscopic factor
    #                 # print(fresco_cross_section_list[2][i])
    #
    #                 dummy_variable = fresco_cross_section_list[0][i]
    #                 fresco_cross_section_list[0][i] = dummy_variable*current_spec/spec_list[0]
    #                 fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]
    #                 if (add_compound_boolean):
    #                     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_sum_list[i]
    #                 # print(fresco_cross_section_list[2][i])
    #                 # print("\n")
    #
    #     #
    #     #
    #     #         print(spec_list[0])
    #
    #
    #
    #
    #     #Fitting L=2 data to points around 30 degrees
    #     if(fit_boolean):
    #         if (spec_list[2]>0.001 and E_x != '7261'):
    #
    #
    #             #Set new limits for spec factors
    #
    #             upper_spec = 1.5*spec_list[2]
    #             lower_spec = 0.5*spec_list[2]
    #
    #
    #
    #             for j in range(0,5): #loop that fits every time
    #                 print("\n\nlower spec is {}, upper spec is {}".format(lower_spec, upper_spec))
    #
    #                 # upper_spec_list[181] = 3.0*fresco_cross_section_list[2]
    #                 # lower_spec_list[181] = lower_spec*fresco_cross_section_list[2]
    #                 for i in range (0,175):
    #                     # print(fresco_cross_section_list[2][i])
    #                     # print(upper_spec_list[i])
    #                     if(add_compound_boolean):
    #                         upper_spec_list[i] = (upper_spec/spec_list[2])*compound_fresco_cs_list[i]
    #                         lower_spec_list[i] = (lower_spec/spec_list[2])*compound_fresco_cs_list[i]
    #                     else:
    #                         upper_spec_list[i] = (upper_spec/spec_list[2])*fresco_cross_section_list[2][i]
    #                         lower_spec_list[i] = (lower_spec/spec_list[2])*fresco_cross_section_list[2][i]
    #
    #                 #Calculate difference between fresco and exp data points for both spec factor points
    #                 sum_diff = 0
    #                 for i in range(6,len(target_4_com_angles_list)): #contribution from target 4
    #                     rounded_angle = round(target_4_com_angles_list[i])
    #                     # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
    #                     rounded_angle = int(rounded_angle)
    #                     # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
    #                     # if(add_compound_boolean):
    #                         # diff = target_4_cs_list[i] -
    #                     diff = target_4_cs_list[i] - upper_spec_list[rounded_angle]
    #                     diff = abs(diff)
    #                     sum_diff = sum_diff + abs(diff)
    #                     # print("diff is {}, giving sum diff of {}".format(diff, sum_diff))
    #
    #                 for i  in range(0,len(target_5_com_angles_list)):
    #                     rounded_angle = round(target_5_com_angles_list[i])
    #                     # print("target 5 value is {} at angle {}, rounded to {}".format(target_5_cs_list[i], target_5_com_angles_list[i], rounded_angle))
    #                     rounded_angle = int(rounded_angle)
    #                     # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
    #                     diff = target_5_cs_list[i] - upper_spec_list[rounded_angle]
    #                     sum_diff = sum_diff + abs(diff)
    #                     # print("sum diff is {}".format(sum_diff))
    #
    #                 upper_sum_diff = sum_diff
    #
    #                 sum_diff = 0
    #
    #                 for i in range(6,len(target_4_com_angles_list)):
    #                     rounded_angle = round(target_4_com_angles_list[i])
    #                     # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
    #                     rounded_angle = int(rounded_angle)
    #                     # print("fresco value at angle {} is {}".format(rounded_angle, lower_spec_list[rounded_angle]))
    #                     diff = target_4_cs_list[i] - lower_spec_list[rounded_angle]
    #                     sum_diff = sum_diff + abs(diff)
    #                     # print("sum diff is {}".format(sum_diff))
    #
    #                 for i  in range(0,len(target_5_com_angles_list)):
    #                     rounded_angle = round(target_5_com_angles_list[i])
    #                     # print("target 5 value is {} at angle {}, rounded to {}".format(target_5_cs_list[i], target_5_com_angles_list[i], rounded_angle))
    #                     rounded_angle = int(rounded_angle)
    #                     # print("fresco value at angle {} is {}".format(rounded_angle, lower_spec_list[rounded_angle]))
    #                     diff = target_5_cs_list[i] - lower_spec_list[rounded_angle]
    #                     sum_diff = sum_diff + abs(diff)
    #                     # print("sum diff is {}".format(sum_diff))
    #
    #
    #                 lower_sum_diff = sum_diff
    #
    #                 print("lower sum diff is {}, upper sum diff is {}".format(lower_sum_diff,upper_sum_diff))
    #
    #                 #Decide which is closer and reset spec factor limits
    #                 if(lower_sum_diff < upper_sum_diff): #lower gives a better estimate
    #                     print("lower gives a better estimate, with a spec factor of {}".format(lower_spec))
    #                     lower_spec = lower_spec
    #                     # upper_spec = (upper_spec + lower_spec)/2.0
    #                     upper_spec = lower_spec + (upper_spec-lower_spec)*0.7
    #                     current_spec = lower_spec
    #
    #                 if(lower_sum_diff > upper_sum_diff): #upper gives a better estimate
    #                     print("upper gives a better estimate, with a sepc factor of {}".format(upper_spec))
    #                     upper_spec = upper_spec
    #                     # lower_spec = (lower_spec + upper_spec)/2.0
    #                     lower_spec = upper_spec - (upper_spec-lower_spec)*0.7
    #                     current_spec = upper_spec
    #
    #             #Outside fitting loop
    #             print("best spectroscopic factor is {}, input spectroscopic factor was {}".format(current_spec, spec_list[2]))
    #             S_2 = current_spec
    #             #changing spectroscopic factor that fresco output is scaled by
    #
    #             for i in range (0,174): #changing model output to new spectroscopic factor
    #                 # print(fresco_cross_section_list[2][i])
    #
    #                 dummy_variable = fresco_cross_section_list[2][i]
    #                 fresco_cross_section_list[2][i] = dummy_variable*current_spec/spec_list[2]
    #                 if (add_compound_boolean):
    #                     if(spec_list[0] < 0.001): #if state is pure l=2 transition, then compound added to just that list
    #                         compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[2][i]
    #                     else:
    #                         fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i];
    #                         compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_sum_list[i]
    #                 # print(fresco_cross_section_list[2][i])
    #                 # print("\n")
    #
    #
    #             print(spec_list[2])


        #Fitting compound nuclear contribution to data
        # if(fit_boolean):
        #     original_state_compound_factor = state_compound_factor
        #     if (add_compound_boolean == 1 and state_compound_factor != 0 and E_x != '7261'):
        #         print("fitting compound contribution")
        #
        #         #Set new limits for spec factors
        #
        #         upper_spec = 1.5*state_compound_factor
        #         lower_spec = 0.5*state_compound_factor
        #
        #
        #
        #         for j in range(0,5): #loop that fits every time
        #             print("\n\nlower state_compound_factor is {}, upper state_compound_factor is {}".format(lower_spec, upper_spec))
        #
        #             for i in range (0,175):
        #                 # print(fresco_cross_section_list[2][i])
        #                 # print(upper_spec_list[i])
        #                 #Making upper and lower angular distributions
        #                 if(add_compound_boolean):
        #                     upper_spec_list[i] = (upper_spec/state_compound_factor)*compound_cs_list[i] + fresco_cross_section_sum_list[i]
        #                     lower_spec_list[i] = (lower_spec/state_compound_factor)*compound_cs_list[i] + fresco_cross_section_sum_list[i]
        #                 else:
        #                     upper_spec_list[i] = (upper_spec/state_compound_factor)*compound_cs_list[i] + fresco_cross_section_list[2][i]
        #                     lower_spec_list[i] = (lower_spec/state_compound_factor)*compound_cs_list[i] + fresco_cross_section_list[2][i]
        #
        #             #Calculate difference between fresco and exp data points for both spec factor points
        #             sum_diff = 0
        #             for i in range(0,len(target_4_com_angles_list)): #contribution from target 4
        #                 rounded_angle = round(target_4_com_angles_list[i])
        #                 # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
        #                 rounded_angle = int(rounded_angle)
        #                 # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
        #                 # if(add_compound_boolean):
        #                     # diff = target_4_cs_list[i] -
        #                 diff = target_4_cs_list[i] - upper_spec_list[rounded_angle]
        #                 diff = abs(diff)
        #                 sum_diff = sum_diff + abs(diff)
        #                 # print("diff is {}, giving sum diff of {}".format(diff, sum_diff))
        #         #
        #             for i  in range(0,len(target_5_com_angles_list)):
        #                 rounded_angle = round(target_5_com_angles_list[i])
        #                 # print("target 5 value is {} at angle {}, rounded to {}".format(target_5_cs_list[i], target_5_com_angles_list[i], rounded_angle))
        #                 rounded_angle = int(rounded_angle)
        #                 # print("fresco value at angle {} is {}".format(rounded_angle, upper_spec_list[rounded_angle]))
        #                 diff = target_5_cs_list[i] - upper_spec_list[rounded_angle]
        #                 sum_diff = sum_diff + abs(diff)
        #                 # print("sum diff is {}".format(sum_diff))
        #
        #             upper_sum_diff = sum_diff
        #
        #             sum_diff = 0
        #
        #             for i in range(0,len(target_4_com_angles_list)):
        #                 rounded_angle = round(target_4_com_angles_list[i])
        #                 # print("target 4 value is {} at angle {}, rounded to {}".format(target_4_cs_list[i], target_4_com_angles_list[i], rounded_angle))
        #                 rounded_angle = int(rounded_angle)
        #                 # print("fresco value at angle {} is {}".format(rounded_angle, lower_spec_list[rounded_angle]))
        #                 diff = target_4_cs_list[i] - lower_spec_list[rounded_angle]
        #                 sum_diff = sum_diff + abs(diff)
        #                 # print("sum diff is {}".format(sum_diff))
        #
        #             for i  in range(0,len(target_5_com_angles_list)):
        #                 rounded_angle = round(target_5_com_angles_list[i])
        #                 # print("target 5 value is {} at angle {}, rounded to {}".format(target_5_cs_list[i], target_5_com_angles_list[i], rounded_angle))
        #                 rounded_angle = int(rounded_angle)
        #                 # print("fresco value at angle {} is {}".format(rounded_angle, lower_spec_list[rounded_angle]))
        #                 diff = target_5_cs_list[i] - lower_spec_list[rounded_angle]
        #                 sum_diff = sum_diff + abs(diff)
        #                 # print("sum diff is {}".format(sum_diff))
        #
        #
        #             lower_sum_diff = sum_diff
        #
        #             print("lower sum diff is {}, upper sum diff is {}".format(lower_sum_diff,upper_sum_diff))
        #
        #             #Decide which is closer and reset spec factor limits
        #             if(lower_sum_diff < upper_sum_diff): #lower gives a better estimate
        #                 print("lower gives a better estimate, with a spec factor of {}".format(lower_spec))
        #                 lower_spec = lower_spec
        #                 # upper_spec = (upper_spec + lower_spec)/2.0
        #                 upper_spec = lower_spec + (upper_spec-lower_spec)*0.7
        #                 current_spec = lower_spec
        #
        #             if(lower_sum_diff > upper_sum_diff): #upper gives a better estimate
        #                 print("upper gives a better estimate, with a sepc factor of {}".format(upper_spec))
        #                 upper_spec = upper_spec
        #                 # lower_spec = (lower_spec + upper_spec)/2.0
        #                 lower_spec = upper_spec - (upper_spec-lower_spec)*0.7
        #                 current_spec = upper_spec
        #
        #         #Outside fitting loop
        #         print("best state_compound_factor factor is {}, input state_compound_factor factor was {}".format(current_spec, original_state_compound_factor))
        #         state_compound_factor = current_spec
        #         #changing state_compound_factor that TALYS output is scaled by
        #
        #         for i in range (0,174): #changing model output to new state_compound_factor
        #             # print(fresco_cross_section_list[2][i])
        #
        #             dummy_variable = compound_cs_list[i]
        #             compound_cs_list[i] = dummy_variable*current_spec/original_state_compound_factor
        #             if (add_compound_boolean):
        #                 if(spec_list[0] < 0.001): #if state is pure l=2 transition, then compound added to just that list
        #                     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[2][i]
        #                 else:
        #                     fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i];
        #                     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_sum_list[i]
        #             # print(fresco_cross_section_list[2][i])
        #             # print("\n")
        #         #
        #         #
        #         #
        #         print(state_compound_factor)
        #         print("j count is {}".format(j))
        # print("k count is {}".format(k))


    if(fitting_uncertainty_boolean):
        print("Adding fitting uncertainty (need to add in 0+ state fitting uncertainties)")
        if(E_x == '5.691'):
            print("For 5.691 MeV state, fitting uncertainty is about 5%")
            for i in range (0, len(target_4_cs_list)):
                curren_uncert = target_4_cs_uncert_list[i]
                print("previous uncertainty is {}".format(curren_uncert))
                new_uncert = curren_uncert*curren_uncert + target_4_cs_list[i]*0.05*target_4_cs_list[i]*0.05
                new_uncert = math.sqrt(new_uncert)
                target_4_cs_uncert_list[i] = new_uncert
                print("new uncertainty is {}\n".format(new_uncert))
            for i in range (0, len(target_5_cs_list)):
                curren_uncert = target_5_cs_uncert_list[i]
                print("previous uncertainty is {}".format(curren_uncert))
                new_uncert = curren_uncert*curren_uncert + target_5_cs_list[i]*0.05*target_5_cs_list[i]*0.05
                new_uncert = math.sqrt(new_uncert)
                target_5_cs_uncert_list[i] = new_uncert
                print("new uncertainty is {}\n".format(new_uncert))

        if(E_x == '6.256'):
            i = 0
            print("For 6.256 MeV state, uncertainties not being read in from file")
            fitting_uncert_file = open("0+_fitting_uncertainties.txt","r")
            fitting_uncert_input = fitting_uncert_file.readlines()
            # for line in fitting_uncert_input:
            #     print(line)
            #     print(lab_angles_list[i])
            #     token = line.split("\t")
            #     if (len(token) > 2):
            #         dummy_angle = token[0]
            #         print(dummy_angle)
            #         dummy_index = target_4_lab_angles_list.index(dummy_angle)

                # if(token[0] == lab_angles_list[i]):
                #     print("match found")
                # i = i + 1
            for i in range (0, len(target_4_cs_list)):
                curren_uncert = target_4_cs_uncert_list[i]
                if target_4_lab_angles_list[i] == 13:
                    fitting_uncert = 0.00
                if target_4_lab_angles_list[i] == 16:
                    fitting_uncert = 0.02
                if target_4_lab_angles_list[i] == 19:
                    fitting_uncert = 0.08
                if target_4_lab_angles_list[i] == 28:
                    fitting_uncert = 0.12
                if target_4_lab_angles_list[i] == 30:
                    fitting_uncert = 0.03
                print("previous uncertainty is {}".format(curren_uncert))
                new_uncert = curren_uncert*curren_uncert + target_4_cs_list[i]*fitting_uncert*target_4_cs_list[i]*fitting_uncert
                new_uncert = math.sqrt(new_uncert)
                target_4_cs_uncert_list[i] = new_uncert
                print("new uncertainty is {}\n".format(new_uncert))
            for i in range (0, len(target_5_cs_list)):
                curren_uncert = target_5_cs_uncert_list[i]
                if target_5_lab_angles_list[i] == 30:
                    fitting_uncert = 0.03
                print("previous uncertainty is {}".format(curren_uncert))
                new_uncert = curren_uncert*curren_uncert + target_5_cs_list[i]*fitting_uncert*target_5_cs_list[i]*fitting_uncert
                new_uncert = math.sqrt(new_uncert)
                target_5_cs_uncert_list[i] = new_uncert
                print("new uncertainty is {}\n".format(new_uncert))

    #merging target 4 and target 5 lists to make fitting simpler

    for i in range (0, len(target_4_cs_list)):
        targets_com_angles_list.append(target_4_com_angles_list[i])
        targets_cs_list.append(target_4_cs_list[i])
        targets_cs_uncertainty_list.append(target_4_cs_uncert_list[i])
    for i in range (0,len(target_5_cs_list)):
        targets_com_angles_list.append(target_5_com_angles_list[i])
        targets_cs_list.append(target_5_cs_list[i])
        targets_cs_uncertainty_list.append(target_5_cs_uncert_list[i])

    print(targets_com_angles_list)
    print(targets_cs_list)
    print(targets_cs_uncertainty_list)

    if(fit_boolean):
        #New fitting method of moving through all available parameter space for L=0, L=2 and compound scaling
        #Setting up flags for different types of states

        #for positive parity states
        if(positive_parity_boolean == 1):
            L_0_flag = L_2_flag = 0
            if(spec_list[0] > 0.001):
                print("L=0 flag set")
                L_0_flag = 1
            if(spec_list[0] < 0.001 and spec_list[2] > 0.001):
                print("just L=2 flag set")
                L_2_flag = 1


            #Get upper and lower limits of compound to iterate through
            min_state_compound_factor = 0.2*state_compound_factor
            max_state_compound_factor = 3.5*state_compound_factor
            print(state_compound_factor, min_state_compound_factor, max_state_compound_factor)

            lower_L_2_spec = 0.5*spec_list[2]
            upper_L_2_spec = 3.5*spec_list[2]
            print(spec_list[2], lower_L_2_spec, upper_L_2_spec)
            if (L_0_flag == 1):
                lower_L_0_spec = 0.2*spec_list[0]
                upper_L_0_spec = 2.5*spec_list[0]
                print(spec_list[0], lower_L_0_spec, upper_L_0_spec)

            #loop that goes through them all and calculates a  difference
            smallest_sum_diff = 5000000000
            best_state_compound_factor = state_compound_factor
            best_L_2_spec = spec_list[2]
            best_L_0_spec = spec_list[0]
            no_fits = 20
            no_L_0_fits = 20


            for i in range (0,no_fits): #loop iterating over every state factor compound
                current_state_compound_factor = min_state_compound_factor + (float(i)/no_fits)*(max_state_compound_factor-min_state_compound_factor)
                # print(current_state_compound_factor)
                for k in range (0, no_fits): #loop iterating over L=2 component
                    current_L_2_spec = lower_L_2_spec + (float(k)/no_fits)*(upper_L_2_spec - lower_L_2_spec)
                    print("current L=2 value is {}".format(current_L_2_spec))
                    for l in range (0, no_L_0_fits): #loop iterating over L=0 component
                        # print(l)
                        if (L_0_flag == 1):
                            current_L_0_spec = lower_L_0_spec + (float(l)/no_L_0_fits)*(upper_L_0_spec - lower_L_0_spec)
                            print("current L=0 value is {}".format(current_L_0_spec))
                            # print(current_L_0_spec)
                        sum_diff = 0
                        for j in range(0,len(targets_com_angles_list)): #loop iterating over every data point to get minimisation
                            rounded_angle = round(targets_com_angles_list[j])
                            # print("target value is {} at angle {}, rounded to {}".format(targets_cs_list[j], targets_com_angles_list[j], rounded_angle))
                            rounded_angle = int(rounded_angle)
                            if(E_x == '3.942' and rounded_angle == 19):
                                continue
                            if(E_x == '5.292' and rounded_angle == 26):
                                continue
                            if(E_x == '5.716' and rounded_angle == 23):
                                continue
                            if(E_x == '6.623' and rounded_angle == 29):
                                continue
                            # print(compound_fresco_cs_list[rounded_angle])
                            # dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + fresco_cross_section_sum_list[rounded_angle]
                            dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + (current_L_2_spec/spec_list[2])*fresco_cross_section_list[2][rounded_angle]
                            if(L_0_flag == 1):
                                dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + (current_L_2_spec/spec_list[2])*fresco_cross_section_list[2][rounded_angle] + (current_L_0_spec/spec_list[0])*fresco_cross_section_list[0][rounded_angle]
                            # print("fresco value at angle {} is {}".format(rounded_angle, dummy_variable))
                            diff = (targets_cs_list[j] - dummy_variable)*(targets_cs_list[j] - dummy_variable)/ (targets_cs_uncertainty_list[j]*targets_cs_uncertainty_list[j])
                            sum_diff = sum_diff + abs(diff)
                            if(j == len(targets_com_angles_list)):
                                print("sum diff is {}\n".format(sum_diff))
                        if(sum_diff < smallest_sum_diff):
                            # print("previous sum_diff was {}, with a state_compound_factor of {} and spec factor of {}".format(smallest_sum_diff, best_state_compound_factor, spec_list[2]))
                            smallest_sum_diff = sum_diff
                            best_state_compound_factor = current_state_compound_factor
                            best_L_2_spec = current_L_2_spec
                            if (L_0_flag == 1):
                                best_L_0_spec = current_L_0_spec
                                print("new sum_diff is {}, with a state_compound_factor of {}, a L=0 spec factor of {}, a L=2 spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_0_spec, best_L_2_spec))
                            # else:
                                # print("new sum_diff is {}, with a state_compound_factor of {} and a spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_2_spec))
                            S_2 = best_L_2_spec
                            if (L_0_flag == 1):
                                S_0 = best_L_0_spec


        #for negative parity states
        if(negative_parity_boolean == 1):
            print("Negative parity boolean turned on")
            L_1_flag = L_3_flag = 0
            if(spec_list[1] > 0.001):
                print("L=1 flag set")
                L_1_flag = 1
            if(spec_list[1] < 0.001 and spec_list[3] > 0.001):
                print("just L=3 flag set")
                L_3_flag = 1


            #Get upper and lower limits of compound to iterate through
            min_state_compound_factor = 0.2*state_compound_factor
            max_state_compound_factor = 3.5*state_compound_factor
            print(state_compound_factor, min_state_compound_factor, max_state_compound_factor)

            lower_L_3_spec = 0.5*spec_list[3]
            upper_L_3_spec = 2.5*spec_list[3]
            print(spec_list[3], lower_L_3_spec, upper_L_3_spec)
            if (L_1_flag == 1):
                lower_L_1_spec = 0.5*spec_list[1]
                upper_L_1_spec = 2.5*spec_list[1]
                print(spec_list[1], lower_L_1_spec, upper_L_1_spec)

            #loop that goes through them all and calculates a  difference
            smallest_sum_diff = 5000000000
            best_state_compound_factor = state_compound_factor
            best_L_3_spec = spec_list[3]
            best_L_1_spec = spec_list[1]
            no_fits = 30
            no_L_1_fits = 30

            for i in range (0,no_fits): #loop iterating over every state factor compound
                current_state_compound_factor = min_state_compound_factor + (float(i)/no_fits)*(max_state_compound_factor-min_state_compound_factor)
                # print(current_state_compound_factor)
                for k in range (0, no_fits): #loop iterating over L=3 component
                    current_L_3_spec = lower_L_3_spec + (float(k)/no_fits)*(upper_L_3_spec - lower_L_3_spec)
                    # print("current L=3 value is {}".format(current_L_3_spec))
                    for l in range (0, no_L_1_fits): #loop iterating over L=1 component
                        # print(l)
                        if (L_1_flag == 1):
                            current_L_1_spec = lower_L_1_spec + (float(l)/no_L_1_fits)*(upper_L_1_spec - lower_L_1_spec)
                            # print("current L=1 value is {}".format(current_L_1_spec))
                            # print(current_L_1_spec)
                        sum_diff = 0
                        for j in range(0,len(targets_com_angles_list)): #loop iterating over every data point to get minimisation
                            rounded_angle = round(targets_com_angles_list[j])
                            # print("target value is {} at angle {}, rounded to {}".format(targets_cs_list[j], targets_com_angles_list[j], rounded_angle))
                            rounded_angle = int(rounded_angle)
                            if(E_x == '3.942' and rounded_angle == 19):
                                continue
                            if(E_x == '5.292' and rounded_angle == 26):
                                continue
                            if(E_x == '5.716' and rounded_angle == 23):
                                continue
                            if(E_x == '6.623' and rounded_angle == 29):
                                continue
                            # print(compound_fresco_cs_list[rounded_angle])
                            # dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + fresco_cross_section_sum_list[rounded_angle]
                            dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + (current_L_3_spec/spec_list[3])*fresco_cross_section_list[3][rounded_angle]
                            if(L_1_flag == 1):
                                dummy_variable = (current_state_compound_factor/state_compound_factor)*compound_cs_list[rounded_angle] + (current_L_3_spec/spec_list[3])*fresco_cross_section_list[3][rounded_angle] + (current_L_1_spec/spec_list[1])*fresco_cross_section_list[1][rounded_angle]
                            # print("fresco value at angle {} is {}".format(rounded_angle, dummy_variable))
                            diff = (targets_cs_list[j] - dummy_variable)*(targets_cs_list[j] - dummy_variable)/ (targets_cs_uncertainty_list[j]*targets_cs_uncertainty_list[j])
                            sum_diff = sum_diff + abs(diff)
                            # print("sum diff is {}\n".format(sum_diff))
                        if(sum_diff < smallest_sum_diff):
                            # print("previous sum_diff was {}, with a state_compound_factor of {} and spec factor of {}".format(smallest_sum_diff, best_state_compound_factor, spec_list[3]))
                            smallest_sum_diff = sum_diff
                            best_state_compound_factor = current_state_compound_factor
                            best_L_3_spec = current_L_3_spec
                            if (L_1_flag == 1):
                                best_L_1_spec = current_L_1_spec
                                print("new sum_diff is {}, with a state_compound_factor of {}, a L=1 spec factor of {}, a L=3 spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_1_spec, best_L_3_spec))
                            # else:
                                # print("new sum_diff is {}, with a state_compound_factor of {} and a spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_3_spec))
                            S_3 = best_L_3_spec
                            if (L_1_flag == 1):
                                S_1 = best_L_1_spec



        if(positive_parity_boolean):
            print("final sum_diff is {}, with a state_compound_factor of {} and a L_2_spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_2_spec))
            if (L_0_flag == 1):
                print("final sum_diff is {}, with a state_compound_factor of {}, L_0_spec factor of {} and a L_2_spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_0_spec, best_L_2_spec))
                #condition to remove compound nuclear component if negligible
                if(L_0_flag == 1):
                    if(compound_cs_list[13]/fresco_cross_section_list[0][13] < 0.01):
                        best_state_compound_factor=0

        if(negative_parity_boolean):
            print("final sum_diff is {}, with a state_compound_factor of {} and a L_3_spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_3_spec))
            if (L_1_flag == 1):
                print("final sum_diff is {}, with a state_compound_factor of {}, L_1_spec factor of {} and a L_3_spec factor of {} ".format(smallest_sum_diff, best_state_compound_factor, best_L_1_spec, best_L_3_spec))
                #condition to remove compound nuclear component if negligible
                if(L_1_flag == 1):
                    if(compound_cs_list[13] < 0.01):
                    # /fresco_cross_section_list[1][13] < 0.01):
                        best_state_compound_factor=0



        if(positive_parity_boolean):
            for i in range (0,170):
                fresco_cross_section_list[2][i] = (best_L_2_spec/spec_list[2])*fresco_cross_section_list[2][i]
                compound_cs_list[i] = (best_state_compound_factor/state_compound_factor)*compound_cs_list[i]
                compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[2][i]
                if(L_0_flag == 1):
                    fresco_cross_section_list[0][i] = (best_L_0_spec/spec_list[0])*fresco_cross_section_list[0][i]
                    compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]
                    fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]

            if(L_0_flag == 1):
                print(fresco_cross_section_list[0][32], fresco_cross_section_list[2][32], compound_cs_list[32], compound_fresco_cs_list[32])
                # compound_fresco_cs_list[i] = (best_state_compound_factor/state_compound_factor)*compound_cs_list[i] + fresco_cross_section_sum_list[i]
                #output whichever one gave the smallest difference


        #Finding differences between experimental results and fresco for L = 2 transitions to

    if(diff_boolean):
        if(L_two > 0.0001):

            exp_fresco_diff_output.write("%s\n" % (E_x))

            for i in range (3,7):
                rounded_angle = round(target_5_com_angles_list[i])
                # print("target 5 value is {} at angle {}, rounded to {}".format(target_5_cs_list[i], target_5_com_angles_list[i], rounded_angle))
                rounded_angle = int(rounded_angle)
                # print("fresco value at angle {} is {}".format(rounded_angle,fresco_cross_section_list[2][rounded_angle]))
                diff = target_5_cs_list[i] - fresco_cross_section_list[2][rounded_angle]
                # print("difference is {}".format(diff))
                exp_fresco_diff_output.write("%f\n" % diff)

    if(OMP_uncertainty_boolean):
        print("adding OMP uncertainty thickness")

        if(L_zero + L_two > 1.1):
            print("OMP uncertainties of l=0 transfer \n")
            #Opening file for l=2 transfer std deviation
            OMP_uncertainty_L_2_input_name = '/home/s1771161/Mg25_d_p_Mg26/src/OMP_stddev_L_2.txt'
            OMP_uncertainty_L_2_input_file = open(OMP_uncertainty_L_2_input_name)
            OMP_uncertainty_L_2_input = OMP_uncertainty_L_2_input_file.readlines()
            for line in OMP_uncertainty_L_2_input: #loop extracts standard deviations for l=2 transfer from file
                token = line.split(" ")
                token[-1] = token[-1].strip() #removes '\n' from final element
                OMP_angle_L_2_list.append(token[0])
                OMP_stddev_L_2_list.append(token[1])

            OMP_uncertainty_L_0_input_name = '/home/s1771161/Mg25_d_p_Mg26/src/OMP_stddev_L_0.txt'
            OMP_uncertainty_L_0_input_file = open(OMP_uncertainty_L_0_input_name)
            OMP_uncertainty_L_0_input = OMP_uncertainty_L_0_input_file.readlines()
            for line in OMP_uncertainty_L_0_input: #loop extracts standard deviations for l=2 transfer from file
                token = line.split(" ")
                token[-1] = token[-1].strip() #removes '\n' from final element
                OMP_angle_L_0_list.append(token[0])
                OMP_stddev_L_0_list.append(token[1])




            for i in range (0,len(targets_com_angles_list)): #loop removing all uncertainties to be added again
                targets_cs_uncertainty_list.pop(0)

            for i in range (0,len(target_4_cs_uncert_list)): #loop adds standard deviations into data
                int_angle = int(target_4_lab_angles_list[i]) #finds angle that matches stddev data
                #Calculate L=0 uncert
                #Calcuate l=2 uncert
                print(fresco_angle_list[0][int_angle], fresco_cross_section_list[0][int_angle], fresco_cross_section_list[2][int_angle], fresco_cross_section_sum_list[int_angle])
                L_0_sum_ratio = fresco_cross_section_list[0][int_angle]/fresco_cross_section_sum_list[int_angle]
                L_2_sum_ratio = fresco_cross_section_list[2][int_angle]/fresco_cross_section_sum_list[int_angle]
                print(L_0_sum_ratio, L_2_sum_ratio, OMP_stddev_L_0_list[int_angle], OMP_stddev_L_2_list[int_angle])

                OMP_uncert = L_0_sum_ratio*float(OMP_stddev_L_0_list[int_angle])*target_4_cs_list[i] + L_2_sum_ratio*float(OMP_stddev_L_2_list[int_angle])*target_4_cs_list[i]

                # OMP_uncert = float(OMP_stddev_L_2_list[int_angle])*target_4_cs_list[i] #scales stddev to current data point
                print("for angle {}, cross section is {}, OMP uncert is {}".format(int_angle,target_4_cs_list[i],OMP_uncert))

                #Get current total uncertainty out of list
                total_uncert = target_4_cs_uncert_list[i]
                print("original total uncert is {}".format(total_uncert))
                # print("total uncert starts as {}".format(total_uncert))
                total_uncert = total_uncert*total_uncert #removes square root
                total_uncert = total_uncert + OMP_uncert*OMP_uncert # adding two uncertainties in quadrature
                total_uncert = math.sqrt(total_uncert) #square rooting again
                # print("adding OMP uncert of {}  gives a total uncert of {}\n".format(OMP_uncert, total_uncert))
                print("new total uncert is {}".format(total_uncert))
                target_4_cs_uncert_list[i] = total_uncert #adding new total uncertainty back into list
                targets_cs_uncertainty_list.append(total_uncert)


            for i in range (0,len(target_5_cs_uncert_list)): #loop adds standard deviations into data
                int_angle = int(target_5_lab_angles_list[i]) #finds angle that matches stddev data
                OMP_uncert = float(OMP_stddev_L_2_list[int_angle])*target_5_cs_list[i] #scales stddev to current data point
                # print("for angle {}, cross section is {}, OMP uncert is {}".format(int_angle,target_5_cs_list[i],OMP_uncert))

                print(fresco_angle_list[0][int_angle], fresco_cross_section_list[0][int_angle], fresco_cross_section_list[2][int_angle], fresco_cross_section_sum_list[int_angle])
                L_0_sum_ratio = fresco_cross_section_list[0][int_angle]/fresco_cross_section_sum_list[int_angle]
                L_2_sum_ratio = fresco_cross_section_list[2][int_angle]/fresco_cross_section_sum_list[int_angle]
                print(L_0_sum_ratio, L_2_sum_ratio)

                OMP_uncert = L_0_sum_ratio*float(OMP_stddev_L_0_list[int_angle])*target_5_cs_list[i] + L_2_sum_ratio*float(OMP_stddev_L_2_list[int_angle])*target_5_cs_list[i]

                #Get current total uncertainty out of list
                total_uncert = target_5_cs_uncert_list[i]
                # print("total uncert starts as {}".format(total_uncert))
                total_uncert = total_uncert*total_uncert #removes square root
                total_uncert = total_uncert + OMP_uncert*OMP_uncert # adding two uncertainties in quadrature
                total_uncert = math.sqrt(total_uncert) #square rooting again
                # print("adding OMP uncert of {}  gives a total uncert of {}\n".format(OMP_uncert, total_uncert))
                target_5_cs_uncert_list[i] = total_uncert #adding new total uncertainty back into list
                targets_cs_uncertainty_list.append(total_uncert)
                #l=2 contribution works, now need to add l=0 contribution

            OMP_uncertainty_L_2_input_file.close()

        else:
            print("Adding OMP uncertainties for just l=2 transfer")
            OMP_uncertainty_input_name = '/home/s1771161/Mg25_d_p_Mg26/src/OMP_stddev_L_2.txt'
            OMP_uncertainty_input_file = open(OMP_uncertainty_input_name)
            OMP_uncertainty_input = OMP_uncertainty_input_file.readlines()
            for line in OMP_uncertainty_input: #loop extracts standard deviations from file
                token = line.split(" ")
                token[-1] = token[-1].strip() #removes '\n' from final element
                OMP_angle_list.append(token[0])
                OMP_stddev_list.append(token[1])
                # print(token)

            #loop calculates ratios for l=0 and l=2 contributions from current state

            # for i in range (0,len(target_4_cs_uncert_list)):  #loop calculates ratios for l=0 and l=2 contributions from current state

            for i in range (0,len(targets_com_angles_list)): #loop removing all uncertainties to be added again
                targets_cs_uncertainty_list.pop(0)

            for i in range (0,len(target_4_cs_uncert_list)): #loop adds standard deviations into data
                int_angle = int(target_4_lab_angles_list[i]) #finds angle that matches stddev data

                OMP_uncert = float(OMP_stddev_list[int_angle])*target_4_cs_list[i] #scales stddev to current data point
                # print("for angle {}, cross section is {}, OMP uncert is {}".format(int_angle,target_4_cs_list[i],OMP_uncert))

                #Get current total uncertainty out of list
                total_uncert = target_4_cs_uncert_list[i]
                # print("total uncert starts as {}".format(total_uncert))
                total_uncert = total_uncert*total_uncert #removes square root
                total_uncert = total_uncert + OMP_uncert*OMP_uncert # adding two uncertainties in quadrature
                total_uncert = math.sqrt(total_uncert) #square rooting again
                # print("adding OMP uncert of {}  gives a total uncert of {}\n".format(OMP_uncert, total_uncert))
                target_4_cs_uncert_list[i] = total_uncert #adding new total uncertainty back into list
                targets_cs_uncertainty_list.append(total_uncert)

            for i in range (0,len(target_5_cs_uncert_list)): #loop adds standard deviations into data
                int_angle = int(target_5_lab_angles_list[i]) #finds angle that matches stddev data
                OMP_uncert = float(OMP_stddev_list[int_angle])*target_5_cs_list[i] #scales stddev to current data point
                # print("for angle {}, cross section is {}, OMP uncert is {}".format(int_angle,target_5_cs_list[i],OMP_uncert))

                #Get current total uncertainty out of list
                total_uncert = target_5_cs_uncert_list[i]
                # print("total uncert starts as {}".format(total_uncert))
                total_uncert = total_uncert*total_uncert #removes square root
                total_uncert = total_uncert + OMP_uncert*OMP_uncert # adding two uncertainties in quadrature
                total_uncert = math.sqrt(total_uncert) #square rooting again
                # print("adding OMP uncert of {}  gives a total uncert of {}\n".format(OMP_uncert, total_uncert))
                target_5_cs_uncert_list[i] = total_uncert #adding new total uncertainty back into list
                targets_cs_uncertainty_list.append(total_uncert)


            OMP_uncertainty_input_file.close()




    #Calculating chi squared
    if(chi_square_boolean):
        print("chi square boolean turned on")
        chi_square = reduced_chi_square = 0
        for i in range(0,len(targets_com_angles_list)):
            if(E_x == '5.716' and i == 3):
                print("\n\n\nSKIPPING\n\n\n")
                continue
            chi_square = chi_square + pow((targets_cs_list[i]-compound_fresco_cs_list[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
            print(targets_cs_list[i], compound_fresco_cs_list[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
        NDF = len(targets_com_angles_list) - 2.0
        fit_reduced_chi_square = chi_square / NDF
        print(NDF)
        chi_square_sd = math.sqrt(2.0*NDF)
        fit_reduced_chi_square_sd = math.sqrt(2.0/NDF)
        print("fit chi square is {} with sd {}, fit reduced chi square is {} with sd {}".format(chi_square, chi_square_sd, fit_reduced_chi_square, fit_reduced_chi_square_sd))
        #Storing chi square calculations

    #Finding values wthin 1sd of chi squared
    if(positive_parity_boolean):

        #Don't change best l=2 specfactors
        if(fit_boolean == 0):
            best_L_2_spec = L_two
        print("Best L=2 specfactor is {}".format(best_L_2_spec))
        if(spec_list[0] > 0.001):
            print("Best L=0 spec factor is {}".format(best_L_0_spec))

        print(len(temp_fresco_cross_section_list[2]))
        for j in range(0,1000): #iterating over different spectrosopic factors
            #varying l=2 contribution
            for i in range (0,170): #iterating through angles
                temp_variable = fresco_cross_section_list[2][i]/best_L_2_spec #brings back to original input - for changing L=2 variable
                # temp_variable = fresco_cross_section_list[0][i]/best_L_0_spec #brings back to original input - for changing L=0 variable
                # print(temp_variable)
                temp_scale = best_L_2_spec*( (float(j)/100) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=2 variable
                # temp_scale = best_L_0_spec*( (float(j)/100) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=0 variable
                temp_fresco_value = temp_scale*temp_variable
                # temp_fresco_cross_section_list[2][i] = best_L_2_spec*(j/100 + 1) * temp_variable
                # compound_cs_list[i] = (best_state_compound_factor/state_compound_factor)*compound_cs_list[i]
                if(spec_list[0] > 0.001):
                    temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[0][i]) #for when changin L=2
                    # temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[2][i]) #for when changin L=0
                else:
                    temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value)
                if(i == 30):
                    print(j, temp_variable, temp_scale, temp_fresco_value)
                    # print(temp_compound_fresco_cs_list[i], compound_cs_list[i], temp_fresco_value)
                # print(temp_compound_fresco_cs_list[i])
                # if(L_0_flag == 1):
                #     fresco_cross_section_list[0][i] = (best_L_0_spec/spec_list[0])*fresco_cross_section_list[0][i]
                #     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]
                #     fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]
            #Calculate chi squared again
            chi_square = reduced_chi_square = 0
            for i in range(0,len(targets_com_angles_list)):
                if(E_x == '5.716' and i == 3 or i== 5):
                    print("\n\n\nSKIPPING\n\n\n")
                    continue
                chi_square = chi_square + pow((targets_cs_list[i]-temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
                print(i,targets_cs_list[i], temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
            NDF = len(targets_com_angles_list) - 2.0
            reduced_chi_square = chi_square / NDF
            # print(NDF)
            chi_square_sd = math.sqrt(2.0*NDF)
            reduced_chi_square_sd = math.sqrt(2.0/NDF)
            print(j)
            print(temp_compound_fresco_cs_list[30])
            print("chi square is {} with sd {}, reduced chi square is {} with sd {}".format(chi_square, chi_square_sd, reduced_chi_square, reduced_chi_square_sd))
            if(reduced_chi_square > (fit_reduced_chi_square + fit_reduced_chi_square_sd)): #When we've gone outside chi square standard deviation
                print("\n\n\n for l=2 component, outside chi squared fit at {}, with a reduced_chi_square of {}, compared to {}  plus {} with {} , giving a diff of {}  \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_2_spec, (temp_scale - best_L_2_spec)))
                # print("\n\n\noutside chi squared fit at {}, with a reduced_chi_square of {}, comapred to {}  plus {} with {} , giving a diff of {} \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_0_spec, (temp_scale - best_L_0_spec)))
                break


    if(negative_parity_boolean):

        #Don't change best L=3 specfactors
        if(fit_boolean == 0):
            best_L_3_spec = L_three
        print("Best L=3 specfactor is {}".format(best_L_3_spec))
        if(spec_list[1] > 0.001):
            print("Best L=1 spec factor is {}".format(best_L_1_spec))

        print(len(temp_fresco_cross_section_list[3]))
        for j in range(0,1000): #iterating over different spectrosopic factors
            #varying L=3 contribution
            for i in range (0,170): #iterating through angles
                # temp_variable = fresco_cross_section_list[3][i]/best_L_3_spec #brings back to original input - for changing L=3 variable
                temp_variable = fresco_cross_section_list[1][i]/best_L_1_spec #brings back to original input - for changing L=1 variable
                # print(temp_variable)
                # temp_scale = best_L_3_spec*( (float(j)/100) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=3 variable
                temp_scale = best_L_1_spec*( (float(j)/100) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=1 variable
                temp_fresco_value = temp_scale*temp_variable
                # temp_fresco_cross_section_list[3][i] = best_L_3_spec*(j/100 + 1) * temp_variable
                # compound_cs_list[i] = (best_state_compound_factor/state_compound_factor)*compound_cs_list[i]
                if(spec_list[1] > 0.001):
                    # temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[1][i]) #for when changin L=3
                    temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[3][i]) #for when changin L=1
                else:
                    temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value)
                if(i == 30):
                    print(j, temp_variable, temp_scale, temp_fresco_value)
                    # print(temp_compound_fresco_cs_list[i], compound_cs_list[i], temp_fresco_value)
                # print(temp_compound_fresco_cs_list[i])
                # if(L_1_flag == 1):
                #     fresco_cross_section_list[1][i] = (best_L_1_spec/spec_list[1])*fresco_cross_section_list[1][i]
                #     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[1][i] + fresco_cross_section_list[3][i]
                #     fresco_cross_section_sum_list[i] = fresco_cross_section_list[1][i] + fresco_cross_section_list[3][i]
            #Calculate chi squared again
            chi_square = reduced_chi_square = 0
            for i in range(0,len(targets_com_angles_list)):
                if(E_x == '5.716' and i == 3 or i== 5):
                    print("\n\n\nSKIPPING\n\n\n")
                    continue
                chi_square = chi_square + pow((targets_cs_list[i]-temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
                print(i,targets_cs_list[i], temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
            NDF = len(targets_com_angles_list) - 2.0
            reduced_chi_square = chi_square / NDF
            # print(NDF)
            chi_square_sd = math.sqrt(2.0*NDF)
            reduced_chi_square_sd = math.sqrt(2.0/NDF)
            print(j)
            print(temp_compound_fresco_cs_list[30])
            print("chi square is {} with sd {}, reduced chi square is {} with sd {}".format(chi_square, chi_square_sd, reduced_chi_square, reduced_chi_square_sd))
            if(reduced_chi_square > (fit_reduced_chi_square + fit_reduced_chi_square_sd)): #When we've gone outside chi square standard deviation
                # print("\n\n\n for L=3 component, outside chi squared fit at {}, with a reduced_chi_square of {}, compared to {}  plus {} with {} , giving a diff of {}  \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_3_spec, (temp_scale - best_L_3_spec)))
                print("\n\n\noutside chi squared fit at {}, with a reduced_chi_square of {}, comapred to {}  plus {} with {} , giving a diff of {} \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_1_spec, (temp_scale - best_L_1_spec)))
                break


    print(reduced_chi_square)
    print('\n\n\n\n\n\n\n\n')
    print(yasue_boolean, E_x)
    # print(best_L_2_spec)
    if(yasue_boolean == 1 and E_x == '5.691'):
        print("Yasue boolean turned on")
        print("Compound factor used is {}".format(best_state_compound_factor))
        for i in range (0,180):
            yasue_com_angle_list.append(i)
            temp_variable = 0.2*fresco_cross_section_list[2][i]/best_L_2_spec #scaling to Yasue's value
            yasue_cs_list.append(temp_variable)
            # print(i,temp_variable)


    #return

    if(both_parities_boolean):
        print("Boolean currently set up to find negative parity limits on 5.7 MeV peak")
        print("Best L=2 spec is {}".format(best_L_2_spec))
        print("Best compound boolean is {}".format(best_state_compound_factor))

        print("Reading in l=1 and l=3 contributions")
        print("Setting up empty arrays")
        l_1_array = []
        l_2_array = []
        l_3_array = []
        compound_array = compound_cs_list

        #Extending compound array to make array 180 elements long to avoid seg faults
        compound_array.append(compound_array[-1])
        compound_array.append(compound_array[-1])
        compound_array.append(compound_array[-1])
        compound_array.append(compound_array[-1])
        total_array = []
        #Reading in l=1 array
        with open("mg26dp_dwba_5.716_3-_1.dat") as l_1_file:
            for counter, line in enumerate(l_1_file):
                # print(line.split("\t")[1])
                l_1_array.append(float(line.split("\t")[1].strip()))
        print("l_1_array:")
        print(l_1_array)

        #Reading in l=2 array
        with open("mg26dp_dwba_5.716_4+_2.dat") as l_2_file:
            for counter, line in enumerate(l_2_file):
                # print(line.split("\t")[1])
                l_2_array.append(float(line.split("\t")[1].strip()))
        print("l_2_array:")
        print(l_2_array)

        #Reading in l=3 array
        with open("mg26dp_dwba_5.716_3-_3.dat") as l_3_file:
            for counter, line in enumerate(l_3_file):
                # print(line.split("\t")[1])
                l_3_array.append(float(line.split("\t")[1].strip()))
        print("l_3_array:")
        print(l_3_array)


        starting_l_1 = 0.001
        starting_l_3 = 0.05
        print("starting_l_1 and starting_l_3 set to {} and {}".format(starting_l_1, starting_l_3))
        # return

        #Finding values wthin 1sd of chi squared
        if(positive_parity_boolean):

            #Don't change best l=2 specfactors

            print("Best L=2 specfactor is {}".format(best_L_2_spec))

            print(len(temp_fresco_cross_section_list[2]))

            #Making fitting to fit l=1, l=2, l=3 (no compound contribution)
            C2S_L1_start, C2S_L1_end = 0.004, 0.018
            C2S_L2_start, C2S_L2_end = 0.03, 0.1
            C2S_L3_start, C2S_L3_end = 0.08, 0.35
            compound_scale_start, compound_scale_end = 0.15, 0.5
            no_C2S = 20

            best_reduced_chi_squared = 1e9

            #Arrays for 3D plotting
            l_1_plotting_array = []
            l_2_plotting_array = []
            l_3_plotting_array = []
            compound_plotting_array = []
            chi_squared_plotting_array = []

            for word in ('testing', 'for', 'loop'):
                print(word)
            print(compound_array)
            print(compound_array[0])
            print(len(compound_array))
            print(len(l_1_array))


            for C2S_L1 in np.arange(C2S_L1_start, C2S_L1_end, (C2S_L1_end - C2S_L1_start)/no_C2S):
                for C2S_L2 in np.arange(C2S_L2_start, C2S_L2_end, (C2S_L2_end - C2S_L2_start)/no_C2S):
                    for C2S_L3 in np.arange(C2S_L3_start, C2S_L3_end, (C2S_L3_end - C2S_L3_start)/no_C2S):
                        for compound_scale in np.arange(compound_scale_start, compound_scale_end, (compound_scale_end - compound_scale_start) / no_C2S):

                            # print(C2S_L1, C2S_L2, C2S_L3, compound_scale)
                            #Calculate total array based on current spec factors
                            # print("Lengths:")
                            # print(len(l_1_array))
                            # print(len(l_2_array))
                            # print(len(l_3_array))
                            # print(len(total_array))
                            # print(l_2_array[0])
                            l_1_plotting_array.append(C2S_L1)
                            l_2_plotting_array.append(C2S_L2)
                            l_3_plotting_array.append(C2S_L3)
                            compound_plotting_array.append(compound_scale)

                            for angle_counter, value in enumerate(l_1_array):
                                # print("Angle counter: {}".format(angle_counter))
                                total_array.append(\
                                C2S_L1 * l_1_array[angle_counter] +\
                                C2S_L2 * l_2_array[angle_counter] +\
                                C2S_L3 * l_3_array[angle_counter] +\
                                compound_scale * compound_array[angle_counter])
                                # print(total_array[angle_counter])

                            #Calculate chi squared again
                            chi_square = reduced_chi_square = 0
                            for i in range(0,len(targets_com_angles_list)):
                                if(E_x == '5.716' and i == 3 or i== 5):
                                    # print("\n\n\nSKIPPING\n\n\n")
                                    continue
                                chi_square = chi_square + pow((targets_cs_list[i]-total_array[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
                                # print(round(targets_com_angles_list[i]), targets_cs_list[i], total_array[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
                                # print(round(targets_com_angles_list[i]), chi_square, targets_cs_list[i], total_array[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i])


                            NDF = len(targets_com_angles_list) - 2.0
                            reduced_chi_square = chi_square / NDF
                            # print(NDF)
                            chi_square_sd = math.sqrt(2.0*NDF)
                            reduced_chi_square_sd = math.sqrt(2.0/NDF)
                            # print("chi square is {} with sd {}, reduced chi square is {} with sd {}".format(chi_square, chi_square_sd, reduced_chi_square, reduced_chi_square_sd))
                            chi_squared_plotting_array.append(reduced_chi_square)
                            if(reduced_chi_square < best_reduced_chi_squared):
                                print("New best reduced chi-squared found: {}\nC2S(l=1): {}\nC2S(l=2): {}\nC2S(l=3): {}\nCompound scale: {}\n".format(reduced_chi_square, C2S_L1, C2S_L2, C2S_L3, compound_scale))
                                best_reduced_chi_squared, best_reduced_chi_square_sd, best_C2S_L1, best_C2S_L2, best_C2S_L3, best_compound_scale = reduced_chi_square, reduced_chi_square_sd, C2S_L1, C2S_L2, C2S_L3, compound_scale
                            del total_array[:]

            print("C2S are being scaled by 3/7 to account for using a 3- state in DWBA calculation, instead of 1-")
            print("Best reduced chi-squared of {} with sd {} found at:\nC2S(l=1): {}\nC2S(l=2): {}\nC2S(l=3): {}\n Compound scale: {}\n".format(best_reduced_chi_squared, best_reduced_chi_square_sd, (7.0 /3.0) * best_C2S_L1, best_C2S_L2, (7.0 /3.0) * best_C2S_L3, best_compound_scale))

            print("Test printing out compound cross section list:")
            print(compound_array)

            print("Making total array to plot:")
            for angle_counter, value in enumerate(l_1_array):
                # print("Angle counter: {}".format(angle_counter))
                total_array.append(\
                best_C2S_L1 * l_1_array[angle_counter] +\
                best_C2S_L2 * l_2_array[angle_counter] +\
                best_C2S_L3 * l_3_array[angle_counter] +\
                best_compound_scale * compound_array[angle_counter])
            chi_square = 0
            for i in range(0,len(targets_com_angles_list)):
                if(E_x == '5.716' and i == 3 or i== 5):
                    # print("\n\n\nSKIPPING\n\n\n")
                    continue
                chi_square = chi_square + pow((targets_cs_list[i]-total_array[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
                # print(round(targets_com_angles_list[i]), targets_cs_list[i], total_array[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
                print(round(targets_com_angles_list[i]), chi_square, targets_cs_list[i], total_array[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i])


            print("Plotting out best chi-squared:")

            plt.plot(best_C2S_L1 * np.array(l_1_array))
            plt.plot(best_C2S_L2 * np.array(l_2_array))
            plt.plot(best_C2S_L3 * np.array(l_3_array))
            plt.plot(best_compound_scale * np.array(compound_array))
            plt.plot(total_array)
            plt.errorbar(target_4_com_angles_list, target_4_cs_list, target_4_cs_uncert_list,label="Target 4",\
            color='firebrick',barsabove=True, marker='.', ms=7, linestyle='None', elinewidth=1, capthick=1,\
            markeredgewidth=2)
            plt.errorbar(target_5_com_angles_list, target_5_cs_list, target_5_cs_uncert_list, label="Target 5",\
            color='firebrick',barsabove=True, marker='.', ms=7, linestyle='None', elinewidth=1, capthick=1, markeredgewidth=2)
            plt.grid(b=True,which='minor',axis='y')
            plt.grid(b=True,which='major',axis='x')
            ax = plt.subplot()
            ax.set_axisbelow(1)
            ax.set_yscale('log')
            ax.tick_params(axis='x',direction = 'in')
            plt.xticks(fontsize=16)
            plt.xlim(0, 60)
            plt.ylim(0.02, 0.5)
            plt.xlabel(r'$\mathrm{\theta_{CoM}}$', fontsize=15)
            plt.ylabel(r'$\mathrm{d\sigma/d\Omega}$ (mb/sr)', fontsize=15,rotation=90, labelpad=-2)
            plt.show()


            fig = plt.figure()
            ax = plt.axes(projection="3d")
            #
            # ax.scatter3D(l_1_plotting_array, l_2_plotting_array, chi_squared_plotting_array,
            # c=l_3_plotting_array, cmap = 'hsv')
            # ax.set_xlabel('L1 C2S')
            # ax.set_ylabel('L2 C2S')
            # ax.set_zlabel('Chi square')
            #
            # plt.show()

            ax.scatter3D(l_1_plotting_array, l_3_plotting_array, chi_squared_plotting_array,
            c=l_2_plotting_array, cmap = 'hsv')
            ax.set_xlabel('L1 C2S')
            ax.set_ylabel('L3 C2S')
            ax.set_zlabel('Chi square')

            plt.show()

            # ax.scatter3D(l_2_plotting_array, l_3_plotting_array, chi_squared_plotting_array,
            # c=l_1_plotting_array, cmap = 'hsv')
            # ax.set_xlabel('L2 C2S')
            # ax.set_ylabel('L3 C2S')
            # ax.set_zlabel('Chi square')

            # plt.show()

            # ax.plot_wireframe(l_1_plotting_array, l_2_plotting_array, chi_squared_plotting_array,
            # c=l_3_plotting_array, cmap = 'hsv')
            #
            # plt.show()


            fig = plt.figure()
            ax = plt.axes(projection="3d")
            # ax.plot_wireframe(np.array(l_1_plotting_array),
            # np.array(l_3_plotting_array), np.array(chi_squared_plotting_array), color='green')
            surf = ax.plot_trisurf(l_1_plotting_array, l_3_plotting_array, chi_squared_plotting_array, linewidth=0, antialiased=False)

            ax.set_xlabel('L=1')
            ax.set_ylabel('L=3')
            ax.set_zlabel('Chi squared')

            plt.show()


            #Need to read in energy difference as a function of angle
            with open("5700keV_energy_diff_per_angle.txt") as energy_diff_file:
                energy_diff_angle_array = []
                energy_diff_array = []
                energy_diff_uncert_array = []

                for line in energy_diff_file:
                    energy_diff_angle, energy_diff, energy_diff_uncert = line.rstrip().split("\t")
                    print(line.rstrip().split("\t"))
                    energy_diff_angle_array.append(float(energy_diff_angle))
                    energy_diff_array.append(float(energy_diff))
                    energy_diff_uncert_array.append(float(energy_diff_uncert))
            print(energy_diff_angle_array)
            print(energy_diff_array)
            print(energy_diff_uncert_array)


            energy_diff_com_angle_array = [com_angles_list[int(angle)-1]  for angle in energy_diff_angle_array]
            plot_energy_diff(energy_diff_angle_array, energy_diff_com_angle_array, energy_diff_array, energy_diff_uncert_array)

            draw_arrows()

            # add_image(fig)
            plt.show()
            # print(com_angles_list)
            # print(energy_diff_angle_array)
            # print(com_angles_list[int(energy_diff_angle_array[0])-1])
            # print(energy_diff_com_angle_array)

            file_string = "5.7MeV_energy_seperation_no_inset.png"

            plt.savefig(file_string,bbox_inches = 'tight',pad_inches = 0)

            return






            return
            for j in range(0,1000): #iterating over different spectrosopic factors
                #varying l=1 contribution
                for i in range (0,170): #iterating through angles
                    temp_variable = l_1_array[i] #brings back to original input
                    # temp_variable = fresco_cross_section_list[0][i]/best_L_0_spec #brings back to original input - for changing L=0 variable
                    # print(temp_variable)
                    temp_scale = starting_l_1*( (float(j)/1000) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=2 variable
                    # temp_scale = best_L_0_spec*( (float(j)/100) + 1) #Change this line to swap between getting higher or lower uncertainties - for changing L=0 variable
                    temp_fresco_value = temp_scale*temp_variable
                    print("Angle: {}, l=1 Xs: {}".format(i, temp_fresco_value))
                    # temp_fresco_cross_section_list[2][i] = best_L_2_spec*(j/100 + 1) * temp_variable
                    # compound_cs_list[i] = (best_state_compound_factor/state_compound_factor)*compound_cs_list[i]
                    if(spec_list[0] > 0.001):
                        temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[0][i]) #for when changin L=1
                        # temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + temp_fresco_value + fresco_cross_section_list[2][i]) #for when changin L=0
                    else:
                        temp_compound_fresco_cs_list[i] = (compound_cs_list[i] + fresco_cross_section_list[2][i] + temp_fresco_value)
                    if(i == 30):
                        print(j, temp_variable, temp_scale, temp_fresco_value)
                        # print(temp_compound_fresco_cs_list[i], compound_cs_list[i], temp_fresco_value)
                    # print(temp_compound_fresco_cs_list[i])
                    # if(L_0_flag == 1):
                    #     fresco_cross_section_list[0][i] = (best_L_0_spec/spec_list[0])*fresco_cross_section_list[0][i]
                    #     compound_fresco_cs_list[i] = compound_cs_list[i] + fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]
                    #     fresco_cross_section_sum_list[i] = fresco_cross_section_list[0][i] + fresco_cross_section_list[2][i]

                #Calculate chi squared again
                chi_square = reduced_chi_square = 0
                for i in range(0,len(targets_com_angles_list)):
                    if(E_x == '5.716' and i == 3 or i== 5):
                        print("\n\n\nSKIPPING\n\n\n")
                        continue
                    chi_square = chi_square + pow((targets_cs_list[i]-temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))]),2)/(targets_cs_uncertainty_list[i]*targets_cs_uncertainty_list[i])
                    print(i,targets_cs_list[i], temp_compound_fresco_cs_list[int(round(targets_com_angles_list[i]))], targets_cs_uncertainty_list[i], chi_square)
                NDF = len(targets_com_angles_list) - 2.0
                reduced_chi_square = chi_square / NDF
                # print(NDF)
                chi_square_sd = math.sqrt(2.0*NDF)
                reduced_chi_square_sd = math.sqrt(2.0/NDF)
                print(j)
                print(temp_compound_fresco_cs_list[30])
                print("chi square is {} with sd {}, reduced chi square is {} with sd {}".format(chi_square, chi_square_sd, reduced_chi_square, reduced_chi_square_sd))
                if(reduced_chi_square > (fit_reduced_chi_square + fit_reduced_chi_square_sd)): #When we've gone outside chi square standard deviation
                    print("\n\n\n for l=1 component, outside chi squared fit at {}, with a reduced_chi_square of {}, compared to {}  plus {} with {} , giving a diff of {}  \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_2_spec, (temp_scale - best_L_2_spec)))
                    # print("\n\n\noutside chi squared fit at {}, with a reduced_chi_square of {}, comapred to {}  plus {} with {} , giving a diff of {} \n\n\n\n".format(temp_scale, reduced_chi_square, fit_reduced_chi_square, fit_reduced_chi_square_sd, best_L_0_spec, (temp_scale - best_L_0_spec)))
                    break









    #Plotting begins
    #Handles defined

    # Target_4, = plt.plot(target_4_lab_angles_list,target_4_cs_list, 'bs', label="Target 4")
    # Target_5, = plt.plot(target_5_lab_angles_list,target_5_cs_list, 'g^', label="Target 5")
    if((L_zero + L_one + L_two + L_three) > 1.1):
        Fresco_Sum, = plt.plot(fresco_angle_sum_list,fresco_cross_section_sum_list,label="Fresco Sum")
    # plt.figure()
    plt.errorbar(target_4_com_angles_list, target_4_cs_list, target_4_cs_uncert_list,label="Target 4",\
    color='firebrick',barsabove=True, marker='.', ms=7, linestyle='None', elinewidth=1, capthick=1,\
    markeredgewidth=2)

    print(len(target_4_cs_list),len(target_4_lab_angles_list),len(target_4_cs_uncert_list))
    plt.errorbar(target_5_com_angles_list, target_5_cs_list, target_5_cs_uncert_list, label="Target 5",\
    color='firebrick',barsabove=True, marker='.', ms=7, linestyle='None', elinewidth=1, capthick=1, markeredgewidth=2)

    if(compound_boolean and E_x != '5.292'):
        Compound_plot, = plt.plot(compound_angles_list, compound_cs_list, ':', markersize = 6,label="Compound")
        Compound_plot.set_markevery=100
        # print(compound_angles_list)
        # print(compound_cs_list)
        print(compound_angles_list[1],compound_cs_list[1])
        print(compound_angles_list[105],compound_cs_list[105])
        # for m in range(0,87):
            # print(compound_angles_list[m],compound_cs_list[m])
            # print(compound_angles_list[m+90],compound_cs_list[m+90])

    if(add_compound_boolean and (yasue_boolean == 0 or E_x != '5.691') and E_x != '5.292'):
        Add_compound_plot, = plt.plot(compound_angles_list, compound_fresco_cs_list, label="DWBA+Compound", color="black")

    if(yasue_boolean == 1 and E_x == '5.691'):
        Yasue, = plt.plot(yasue_com_angle_list, yasue_cs_list, label="Yasue", color="black")
        Yasue.set_dashes([10,2])

    if L_list[0] == 1:
        Fresco_L_0, = plt.plot(fresco_angle_list[0],fresco_cross_section_list[0], label="Fresco L=0", color="mediumpurple")
    if L_list[1] == 1:
        Fresco_L_1, = plt.plot(fresco_angle_list[1],fresco_cross_section_list[1], label="Fresco L=1")
    if L_list[2] == 1:
        if(E_x == '5.292'):
            Fresco_L_2, = plt.plot(fresco_angle_list[2],fresco_cross_section_list[2], label="Fresco L=2", color="black")
        else:
            Fresco_L_2, = plt.plot(fresco_angle_list[2],fresco_cross_section_list[2], label="Fresco L=2", color="darkorange")
    if L_list[3] == 1:
        Fresco_L_3, = plt.plot(fresco_angle_list[3],fresco_cross_section_list[3], label="Fresco L=3")

    # if((L_zero + L_one + L_two + L_three) > 1.1 ):
    if(1==1):
        if L_list[0] == 1:
            Fresco_L_0.set_dashes([2, 2, 10, 2])
        if L_list[1] == 1:
            Fresco_L_1.set_dashes([2, 2, 10, 2])
        if (L_list[2] == 1 and E_x != '5.292'):
            Fresco_L_2.set_dashes([2, 2, 10, 2])
        if L_list[3] == 1:
            Fresco_L_3.set_dashes([2, 2, 10, 2])
    #Set compound contribution to be dashed
    # Compound_plot.set_dashes([2, 2, 10, 2])


    # Fresco, = plt.plot(fresco_angle_list,fresco_cross_section_list, label="Fresco")
    # plt.legend(handles=[Target_4,Target_5, Fresco_Sum])
    plt.xlim(0,60)

    max_value = 3.0*max(target_4_cs_list)
    print(max_value)
    min_value = 0.2*min(target_5_cs_list)
    if(L_list[2] > 0.00001):
        if((min(compound_cs_list) < min_value) and state_compound_factor > 0.001 and compound_cs_list[13]/fresco_cross_section_list[2][13] > 0.01):
            min_value = min(compound_cs_list)
    print(min_value)

    if(yasue_boolean == 1 and E_x == '5.691'):
        max_value = 0.3

    plt.ylim(min_value,max_value)


    # plt.ylabel('cross section [mb/sr]')
    print('E_xis curently {}'.format(E_x))
    if(E_x == '3.588' or E_x == '4.350' or E_x == '4.972'or E_x == '5.691'or E_x == '6.256'or E_x == '6.876' or E_x == '7.100' or E_x == '7.349' or paper_boolean == 1):
        print('E_xis currently {}'.format(E_x))

        # plt.ylabel(r'$\frac{d\sigma}{d\Omega}\left(\frac{mb}{sr}\right)$', fontsize=25,rotation=90, labelpad=-5)
        plt.ylabel(r'$\mathrm{d\sigma/d\Omega}$ (mb/sr)', fontsize=15,rotation=90, labelpad=-2)
    # plt.xlabel(r'angle$\left(^{\circ}\right)$')
    if(paper_boolean == 0 or E_x == '5.691'): #adds theta com axis label to bottom of plot
        plt.xlabel(r'$\mathrm{\theta_{CoM}}$', fontsize=15)

    # if(paper_boolean == 1 and E_x !='6.256'): #removes numbers from x axis
        # x = np.array([0,10,20,30])
        # my_xticks = ['','','','']
        # plt.xticks([x, my_xticks])
        # plt.xticks([])
        # plt.axes.set_axisbelow()
        # ax.xaxis.get_major_ticks.tick.label.set_fontsize(100)


    ax = plt.subplot()
    ax.set_axisbelow(1)
    ax.set_yscale('log')
    ax.tick_params(axis='x',direction = 'in')
    plt.xticks(fontsize=16)
    if(paper_boolean == 1 and E_x !='5.691'): #shrinks numbers from x axis for all but final plot
        plt.xticks(fontsize=0.000000000)
        ax.xaxis.label.set_color('white')

    plt.yticks(fontsize=16)
    # ax.tick_params(axis='both', which='major', labelsize=15)
    if(paper_boolean == 0):
        ax.legend(framealpha=1)
        # markers.legend()
        plt.grid(b=True,which='minor',axis='y')
        plt.grid(b=True,which='major',axis='x')
    title_string = E_x + " MeV " + J_pi + " state"
    print(title_string)

    hfont = {'fontname':'Helvetica'}
    if(paper_boolean == 1): #Adds labels on to plots

        J_pi_string = J_pi[0] + "$^{+}$"
        text_string = J_pi_string + " " + E_x + " MeV"
        if(E_x =='5.292'):
            plt.text(0.6, 0.9, '(b) ' + text_string, fontsize = 15, fontname='Helvetica', horizontalalignment='left',verticalalignment='center',  transform=ax.transAxes)
        if(E_x =='5.691'):
            plt.text(0.5, 0.2, '(d) ' +text_string, fontsize = 15, fontname='Helvetica', horizontalalignment='left',verticalalignment='center',  transform=ax.transAxes)
        if(E_x =='6.125'):
            plt.text(0.6, 0.9, '(a) ' +text_string, fontsize = 15, fontname='Helvetica', horizontalalignment='left',verticalalignment='center',  transform=ax.transAxes)
        if(E_x =='6.256'):
            plt.text(0.6, 0.9, '(c) ' +text_string, fontsize = 15, fontname='Helvetica', horizontalalignment='left',verticalalignment='center',  transform=ax.transAxes)



    if(paper_boolean == 0):
        plt.title(title_string)
    spec_annotation = "$\mathbf {(2J+1)S}$\n""L=0: " + str((2*J+1)*S_0) + "\nL=1: " + str((2*J+1)*S_1) + "\nL=2: " + str((2*J+1)*S_2) + "\nL=3: " + str((2*J+1)*S_3)
    print("S 0 is {}, S 2 is {}".format(S_0, S_2))
    if(paper_boolean == 0):
        plt.text(0.2, 0.2, spec_annotation, horizontalalignment='left',verticalalignment='center', bbox=dict(facecolor='white', alpha=0.9), transform=ax.transAxes)

    fig_size = plt.rcParams["figure.figsize"]
    print("fig siz is {}".format(fig_size))
    #Changing file name because latex doesn't like '.' in middle of file name
    E_x = float(E_x)*1000.0
    E_x = int(E_x)
    E_x = str(E_x)
    file_string = "results_plts/" + E_x + "_keV_" + J_pi + "state.pdf"
    # file_string = "results_plots/" + E_x + "_keV_" + J_pi + "_state.pdf"

    if(compound_boolean):
        file_string = "results_plots_show_compound/" + E_x + "_keV_" + J_pi + "_state.pdf"
    if(add_compound_boolean):
        file_string = "results_plots_bcfec_compound/" + E_x + "_keV_" + J_pi + "_state_combined.pdf"
    if(combined_fit):
        file_string = "results_plots/" + E_x + "_keV_" + J_pi + "_state_combined.pdf"
    if (save_fig):
        plt.savefig(file_string, bbox_inches='tight',pad_inches = 0)
        print("plot saved to {}".format(file_string))
    else:
        plt.show()
        print('else')

    plt.clf()


    print("Fresco values read in from {}".format(fresco_file_name))
    print("Cross-sections read in from {}".format(result_file_name))
    print("\n\n\n")

    # print(fresco_angle_list[2],fresco_cross_section_list[2])
    # for i in range (0, len(target_4_com_angles_list)):
        # print(target_4_com_angles_list[i], target_4_cs_list[i], target_4_cs_uncert_list[i])

    # for i in range (0, len(target_5_com_angles_list)):
        # print(target_5_com_angles_list[i], target_5_cs_list[i], target_5_cs_uncert_list[i])

    # print(compound_angles_list, compound_fresco_cs_list)

    # print(fresco_angle_sum_list)
    # print(fresco_angle_list[1])
    # print(fresco_angle_list[3])

    # print(fresco_cross_section_sum_list[3])
    # print(fresco_cross_section_list[1][3])
    # printfresco_cross_section_list[2][3])

    print(target_4_cs_list)

    return plt

main()

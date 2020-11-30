#Program takes fresco .out output as its input and converts to just angle and angular cross-section


import sys
import matplotlib.pyplot as plt
# import matplotlib.axes.Axes as axes

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

angle_list=[]
# cross_section_list=[]
burlein_angles=[]
burlein_cross_sections=[]
tunl_angles=[]
tunl_cross_sections=[]

def main(): #if running with run_fresco.py, place filename as an argument here

    filename=sys.argv[1]        #if running with run_fresco.py, comment this out #takes in command line argument as input file in " "
    print(filename)

    #Initialising 2-D list
    n=5
    m=180
    cross_section_list=[[] * m for i in range(n)]
    angle_list=[[] * m for i in range(n)]

    #Reading in Burlein (1984) data
    # burlein_in = open("burlein_data.dat","r")
    # burlein_in = open("burlein_data_0+_3.588.dat","r")
   # burlein_in = open("burlein_data_2+_5.289.dat","r")
    # burlein_in = open("burlein_data_4+_5.474.dat","r")

    #burlein_input = burlein_in.readlines()

    #for line in burlein_input:
        # print(line)
     #   token=line.split("\t")
        # print(token)
      #  token[0]=float(token[0])
       # token[3]=float(token[3])
        # print(token[0])
        #burlein_angles.append(token[0])
        #burlein_cross_sections.append(token[3])

    # tunl_in = open("tunl_data_4+_5.474.dat")
  #  tunl_in = open("tunl_data_2+_5.292.dat")
    # tunl_in = open("tunl_data_0+_3.588.dat")

   # tunl_input = tunl_in.readlines()

    #for line in tunl_input:
        # print(line)
     #   token=line.split("\t")
        # print(token)
      #  token[0]=float(token[0])
       # token[1]=float(token[1])
        #tunl_angles.append(token[0])
       # tunl_cross_sections.append(token[1])




    output_filename = filename[:-3] + 'dat'
    print(output_filename)
    filein = open(filename,"r")
    input_data = filein.readlines()    #Read in whole file

    data_flag=0
    reaction_counter=0

    for line in input_data: #reading in data loop begins


        if line.startswith("     2:  J= 0.5+"):
            print("line with info found")
            token=line.split(" ")

            token=filter(None,token)
            #
            print(token)
            print(token[12],token[16])
            J_pi=token[12]
            E_x=token[16]
            E_x=E_x[:6]
            J=J_pi[:1]
            pi=J_pi[-1]
            # token[12]=float(token[12])
            # token[15]=float(token[15])
            # E_x=round_sig(token[15],4)
            print(E_x)
            print(J)
            print(pi)
            # output_filename =
            print("Ex= {} J_pi= {} {}".format(E_x,J,pi))




        if line.startswith("0   3: "): #this finds the third single-particle form factor in the fresco_outputs
                                        #output, which is what it uses as its second output, so this shouldn't need to be changed for different L transfers
            print("line with L transfer number found")
            token=line.split(" ")

            token=filter(None,token)

            print(token)
            print(token[7])
            delta_L=token[7]
            print("L={}".format(delta_L))



            output_filename = "fresco_outputs_bcfec/mg26dp_dwba_" + E_x + "_" + J + pi + "_"+delta_L + ".dat"
            print("output filename is {}".format(output_filename))

    fileout = open(output_filename,"w")

    filein.close()
    filein = open(filename,"r")

    for line in input_data: #reading in data loop begins

        if line.startswith(" Finished all xsecs"):
            data_flag=0
            print("End of all cross-section data reached, data flag set to 0")
            # print(line)



        if data_flag==1:
            token=line.split(" ")
            token=filter(None,token)

            # print(token)
            if hasNumbers(token[0]): #checking if first element actually contains a number, to avoid conversion error

                token[0]=float(token[0])
                token[4]=float(token[4])
                # fileout.write("{}\t{}\n".format(token[0],token[4]))
                angle_list[reaction_counter].append(token[0])
                cross_section_list[reaction_counter].append(token[4])



        if line.startswith (" CROSS SECTIONS FOR OUTGOING p "):
            data_flag=1
            print("Cross-section data reached, data flag set to 1")
            print(line)

        if line.startswith ("   180.00 deg"):
            print("End of one reaction reached")
            print("Reaction counter is {}".format(reaction_counter))
            reaction_counter=reaction_counter+1
            # plt.plot(angle_list,cross_section_list)
            # plt.plot(burlein_angles,burlein_cross_sections,'ro')
            # plt.xlim(0,80)
            # plt.ylabel('cross section [mb/sr]')
            # plt.xlabel('angle/degrees')
            # ax = plt.subplot()
            # ax.set_yscale('log')
            # plt.grid()
            # plt.show()


    #scaling to S-factors from Burlein (1984)
    for i in range (0,m):
        # temp_variable=cross_section_list[i]
        # print(temp_variable)
        cross_section_list[0][i]=1*cross_section_list[0][i]
        cross_section_list[1][i]=1*cross_section_list[1][i]
        print("{}\t{}\n".format(i,cross_section_list[1][i]))
        fileout.write("{}\t{}\n".format(i,cross_section_list[1][i]))

        # cross_section_list[3][i]
        # temp_variable=cross_section_list[0][i] + cross_section_list[2][i]
        # cross_section_list[3].append(temp_variable)


    # print(cross_section_list[0])
    # print(cross_section_list[1])
    print(burlein_cross_sections)
    angle_list[3]=angle_list[1]

    #plotting loop
    # for i in range(0,reaction_counter+1):
    #     print(len(cross_section_list[i]))
    #
    #     plt.plot(angle_list[i],cross_section_list[i])
    #     # print(angle_list[i])
    #     print(cross_section_list[i])
    #

    # plt.plot(angle_list[0],cross_section_list[0])
    #plt.plot(angle_list[1],cross_section_list[1])
    # plt.plot(burlein_angles,burlein_cross_sections,'ro')
    # plt.plot(tunl_angles,tunl_cross_sections,'b+')

    # plt.plot(angle_list[2],cross_section_list[2])
    # plt.plot(angle_list[3],cross_section_list[3])

    #plt.xlim(0,80)
    #plt.ylabel('cross section [mb/sr]')
    #plt.xlabel('angle/degrees')
    #ax = plt.subplot()
    #ax.set_yscale('log')
    #plt.grid()
    #plt.show()

    # print(angle_list[1])
    # print(cross_section_list[1])

    # for i in range (0,180):
    #     fileout.write(i)

    filein.close()         # Good practice
    fileout.close()
    #burlein_in.close()
   # tunl_in.close()
    print(cross_section_list[1][0])
    print(output_filename)
    return()

main() #if running with run_fresco.py, comment this line out

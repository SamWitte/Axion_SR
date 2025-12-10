import numpy as np


def erg_shift_1(n, alph=0.1):
    return 1.0 * (1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4))
#    return 1.0 * (1.0 - alph**2 / (2 * n**2) )

Nmax = 8
file_output = "load_rate_input_Nmax_{:.0f}.txt".format(Nmax)
n_levels = [i for i in range(2,Nmax+1)]
# n_levels = [2, 3, 4, 5]
hold_out_1 = []
hold_out_2 = []
hold_out_3 = []
hold_out_4 = []


# ignore_lvls = ["|321>", "|421>", "|431>", "|432>", "|521>", "|531>", "|541>", "|532>", "|542>", "|543>"]
ignore_lvls = [] # if you want to pre-cut, but i wont do it for now

final_keep = []
final_keep_sve = []
print_all = False
sve_all = True

maxInt = Nmax - 1
# compute all relevant permutations up to nMax
for n1 in n_levels:
    for n2 in n_levels:
        for n3 in n_levels:
            for l1 in range(1, n1):
                for l2 in range(1, n2):
                    for l3 in range(1, n3):
                        
                        for m1 in range(1, l1+1):
                            for m2 in range(1, l2+1):
                                for m3 in range(1, l3+1):
                                    if (m1>l1) or (m2>l2) or (m3>l3) or (l1 >= n1) or (l2 >= n2) or (l3 >= n3):
                                        continue
                                        
                                    # erg
                                    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)
                                    
                                    tag1 = "|{}{}{}>".format(n1,l1,m1)
                                    tag2 = "|{}{}{}>".format(n2,l2,m2)
                                    tag3 = "|{}{}{}>".format(n3,l3,m3)
                                    if (tag1 in ignore_lvls)or(tag2 in ignore_lvls)or(tag3 in ignore_lvls):
                                        continue
                                        
                                    if erg_diff > 1.0: # emission to infinity
                                        # print(n1, n2, l1, l2, erg_diff, test1)
                                        lf = l1 + l2 - l3
                                        mf = m1 + m2 - m3
                                        
                                        tag4 = "Inf ({}{})".format(lf, mf)
                                        
                                        
                                        in_list = False
                                        for i in range(len(hold_out_1)):
                                            if ((hold_out_1[i]==tag1)and(hold_out_2[i]==tag2)and(hold_out_3[i]==tag3) or (hold_out_1[i]==tag2)and(hold_out_2[i]==tag1)and(hold_out_3[i]==tag3)):
                                                in_list=True
                                        
                                        if not in_list:
                                            hold_out_1.append(tag1)
                                            hold_out_2.append(tag2)
                                            hold_out_3.append(tag3)
                                            hold_out_4.append(tag4)
                                            if print_all:
                                                print(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                                            else:
                                                final_keep.append(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                                                outtag = "Inf"
                                                final_keep_sve.append(tag1[1:4] + "    " + tag2[1:4] + "    " + tag3[1:4] + "    " + outtag)

                                        
                                    else: # m sum
                                        test_m = ((m1 + m2 - m3) == 0)
                                        for nfinal in n_levels:
                                            if (tag1 == tag3 or tag2 == tag3) or test_m or nfinal <= 2:
                                                continue
                                            if (erg_diff - erg_shift_1(nfinal)) == 0:
                                                print("Resonant. \t", tag1, " ", tag2, " ", tag3, " nfinal", nfinal)
                                        
                                        
                                        test_l = ((l1 + l2 - l3) == 0)
                                        # l even
                                        test1 = ((l1 + l2 + l3) % 2 ==0)
                                        if test1 and test_m and test_l:
                                            tag4 = "BH"
                                            in_list = False
                                            for i in range(len(hold_out_1)):
                                                if ((hold_out_1[i]==tag1)and(hold_out_2[i]==tag2)and(hold_out_3[i]==tag3) or (hold_out_1[i]==tag2)and(hold_out_2[i]==tag1)and(hold_out_3[i]==tag3)):
                                                    in_list=True
                                            
                                            if not in_list:
                                                hold_out_1.append(tag1)
                                                hold_out_2.append(tag2)
                                                hold_out_3.append(tag3)
                                                hold_out_4.append(tag4)
                                                if print_all:
                                                    print(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                                                else:
                                                    final_keep.append(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                                                    
                                                    outtag = "BH"
                                                
                                                    final_keep_sve.append(tag1[1:4] + "    " + tag2[1:4] + "    " + tag3[1:4] + "    " + outtag)

### here is how we truncate n-expansion
# Store the extra states for inverse processes
extra_states = []  # Will store (n1, l1, m1, n3, l3, m3) tuples

for n1 in n_levels:
    for l1 in range(0, n1):
        for m1 in range(1, l1 + 1):
            if (m1>l1) or (m2>l2) or (m3>l3) or (l1 >= n1) or (l2 >= n2) or (l3 >= n3):
                continue

            n2 = n1
            l2 = l1
            m2 = m1

            m3 = m1 + m2
            l3 = m3
            n3 = l3 + 1

            # Store for inverse process (will create n3 x n3 -> n1 x Inf)
            extra_states.append((n1, l1, m1, n3, l3, m3))

            # erg
            erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)
            
            
            tag1 = "|{}{}{}>".format(n1,l1,m1)
            tag2 = "|{}{}{}>".format(n2,l2,m2)
            tag3 = "|{}{}{}>".format(n3,l3,m3)
            if (tag1 in ignore_lvls)or(tag2 in ignore_lvls)or(tag3 in ignore_lvls):
                continue
            
            
            if erg_diff > 1.0: # emission to infinity
                # print(n1, n2, l1, l2, erg_diff, test1)
                lf = l1 + l2 - l3
                mf = m1 + m2 - m3
                
                tag4 = "Inf ({}{})".format(lf, mf)
                
                
                in_list = False
                for i in range(len(hold_out_1)):
                    if ((hold_out_1[i]==tag1)and(hold_out_2[i]==tag2)and(hold_out_3[i]==tag3) or (hold_out_1[i]==tag2)and(hold_out_2[i]==tag1)and(hold_out_3[i]==tag3)):
                        in_list=True
                
                if not in_list:
                    hold_out_1.append(tag1)
                    hold_out_2.append(tag2)
                    hold_out_3.append(tag3)
                    hold_out_4.append(tag4)
                    if print_all:
                        print(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                    else:
                        final_keep.append(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                        outtag = "Inf"
                        final_keep_sve.append(tag1[1:4] + "    " + tag2[1:4] + "    " + tag3[1:4] + "    " + outtag)

                
            else: # m sum
                test_m = ((m1 + m2 - m3) == 0)
                test_l = ((l1 + l2 - l3) == 0)
                # l even
                test1 = ((l1 + l2 + l3) % 2 ==0)
                if test1 and test_m and test_l:

                    tag4 = "BH"
                    in_list = False
                    for i in range(len(hold_out_1)):
                        if ((hold_out_1[i]==tag1)and(hold_out_2[i]==tag2)and(hold_out_3[i]==tag3) or (hold_out_1[i]==tag2)and(hold_out_2[i]==tag1)and(hold_out_3[i]==tag3)):
                            in_list=True
                    
                    if not in_list:
                        hold_out_1.append(tag1)
                        hold_out_2.append(tag2)
                        hold_out_3.append(tag3)
                        hold_out_4.append(tag4)
                        if print_all:
                            print(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                        else:
                            final_keep.append(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                            
                            outtag = "BH"
                        
                            final_keep_sve.append(tag1[1:4] + "    " + tag2[1:4] + "    " + tag3[1:4] + "    " + outtag)


### Add inverse processes for extra states beyond Nmax
# If we have n_low x n_low -> n_high x BH, add n_high x n_high -> n_low x Inf
for (n_low, l_low, m_low, n_high, l_high, m_high) in extra_states:
    # Create inverse process: n_high x n_high -> n_low x Inf
    n1 = n_high
    l1 = l_high
    m1 = m_high

    n2 = n_high
    l2 = l_high
    m2 = m_high

    n3 = n_low
    l3 = l_low
    m3 = m_low

    # Check energy - should emit to infinity
    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)

    tag1 = "|{}{}{}>".format(n1,l1,m1)
    tag2 = "|{}{}{}>".format(n2,l2,m2)
    tag3 = "|{}{}{}>".format(n3,l3,m3)

    if erg_diff > 1.0:  # emission to infinity
        lf = l1 + l2 - l3
        mf = m1 + m2 - m3

        tag4 = "Inf ({}{})".format(lf, mf)

        # Check if already in list
        in_list = False
        for i in range(len(hold_out_1)):
            if ((hold_out_1[i]==tag1)and(hold_out_2[i]==tag2)and(hold_out_3[i]==tag3) or (hold_out_1[i]==tag2)and(hold_out_2[i]==tag1)and(hold_out_3[i]==tag3)):
                in_list=True

        if not in_list:
            hold_out_1.append(tag1)
            hold_out_2.append(tag2)
            hold_out_3.append(tag3)
            hold_out_4.append(tag4)
            if print_all:
                print(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
            else:
                final_keep.append(tag1 + " x " + tag2 + " ---> " + tag3 + " x " + tag4)
                outtag = "Inf"
                final_keep_sve.append(tag1[1:4] + "    " + tag2[1:4] + "    " + tag3[1:4] + "    " + outtag)


if sve_all:
    with open(file_output, 'w') as f:
        for i in range(len(hold_out_1)):
            f.write(final_keep_sve[i] + " \n")  # python will convert \n to os.linesep
    
            

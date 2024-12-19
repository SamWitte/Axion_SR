import numpy as np


def erg_shift_1(n, alph=0.1):
    return 1.0 * (1.0 - alph**2 / (2 * n**2) - alph**4 / (8 * n**4))

n_levels = [2, 3, 4, 5]
# n_levels = [2, 3, 4, 5]
hold_out_1 = []
hold_out_2 = []
hold_out_3 = []
hold_out_4 = []

# [311 cannot grow bc 211 x 311 -> 322 BH]
# ignore_lvls = ["|321>", "|421>", "|431>", "|432>"]
ignore_lvls = ["|321>", "|421>", "|431>", "|432>", "|521>", "|531>", "|541>", "|532>", "|542>", "|543>"]
# ignore_lvls = []
keep_lvls = ["|511>"]
# keep_lvls = []
final_keep = []
final_keep_sve = []
print_all = False
sve_all = True

maxInt = np.max(n_levels) - 1

for n1 in n_levels:
    for n2 in n_levels:
        for n3 in n_levels:
            for l1 in range(0, maxInt + 1):
                for l2 in range(0, maxInt + 1):
                    for l3 in range(0, maxInt + 1):
                        
                        for m1 in range(1, maxInt + 1):
                            for m2 in range(1, maxInt + 1):
                                for m3 in range(1, maxInt + 1):
                                    if (m1>l1) or (m2>l2) or (m3>l3) or (l1 >= n1) or (l2 >= n2) or (l3 >= n3):
                                        continue
                                        
                                  
                                    
                                    # erg
                                    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)
                                    print(erg_diff)
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


if len(keep_lvls) > 0:
    
        if (keep_lvls[0] == hold_out_1[i]) or (keep_lvls[0] == hold_out_2[i]) or (keep_lvls[0] == hold_out_3[i]):
            print(final_keep[i])
if sve_all:
    with open('load_rate_input.txt', 'w') as f:
        for i in range(len(hold_out_1)):
            f.write(final_keep_sve[i] + " \n")  # python will convert \n to os.linesep
    
            

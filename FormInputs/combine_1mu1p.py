import os

josh_list  = [ "FinalSelection0_overlay_run1.txt", "FinalSelection0_intrinsics_run1.txt", "FinalSelection0_ncpi0_run1.txt", "FinalSelection0_ccpi0_run1.txt", "FinalSelection0_ext_run1.txt",
               "FinalSelection0_overlay_run2.txt", "FinalSelection0_intrinsics_run2.txt", 
               "FinalSelection0_overlay_run3.txt", "FinalSelection0_intrinsics_run3.txt", "FinalSelection0_ncpi0_run3.txt", "FinalSelection0_ccpi0_run3.txt", "FinalSelection0_ext_run3.txt" ]
katie_list = [ "May05_pi0_weights_1mu1p_numu_run1.txt", "May05_pi0_weights_1mu1p_nue_run1.txt", "May05_pi0_weights_1mu1p_ncpi0_run1.txt", "May05_pi0_weights_1mu1p_ccpi0_run1.txt", "empty.txt",
               "May05_pi0_weights_1mu1p_numu_run2.txt", "May05_pi0_weights_1mu1p_nue_run2.txt", 
               "May05_pi0_weights_1mu1p_numu_run3.txt", "May05_pi0_weights_1mu1p_nue_run3.txt", "May05_pi0_weights_1mu1p_ncpi0_run3.txt", "May05_pi0_weights_1mu1p_ccpi0_run3.txt", "empty.txt" ]
out_list   = [ "sel1mu1p_bnb_run1_eventlist.txt", "sel1mu1p_nue_run1_eventlist.txt", "sel1mu1p_ncpi0_run1_eventlist.txt", "sel1mu1p_ccpi0_run1_eventlist.txt", "sel1mu1p_extbnb_run1_eventlist.txt",
               "sel1mu1p_bnb_run2_eventlist.txt", "sel1mu1p_nue_run2_eventlist.txt",
               "sel1mu1p_bnb_run3_eventlist.txt", "sel1mu1p_nue_run3_eventlist.txt", "sel1mu1p_ncpi0_run3_eventlist.txt", "sel1mu1p_ccpi0_run3_eventlist.txt", "sel1mu1p_extbnb_run3_eventlist.txt" ]

top_dir   = "/uboone/app/users/yatesla/othersys_mcc9/form_input_to_sbnfit"
josh_dir  = os.path.join(top_dir, "fromJosh")
katie_dir = os.path.join(top_dir, "fromKatie")

for i in range(len(out_list)):
        
    print "for file %s..." % out_list[i]
    
    # Read in the pi0 weights
    pi0_weight_dict = {}
    with open( os.path.join(katie_dir, katie_list[i]), 'r' ) as k:
        for l in k:
            if l.split(',')[0] == "run": continue
            cols = l.strip().split(',')
            run = int(cols[0])
            subrun = int(cols[1])
            event = int(cols[2])
            energy = float(cols[3])
            weight = float(cols[4])
            rse = tuple( (run, subrun, event, energy) )
            pi0_weight_dict[rse] = weight
        
    # Write out the file with the reco vars plus the weight
    tot_ctr = 0
    wgh_ctr = 0
    err_ctr = 0
    with open( os.path.join(top_dir, out_list[i]), 'w' ) as out:
        with open( os.path.join(josh_dir, josh_list[i]), 'r' ) as j:
            for l in j:
                if l.startswith('run'): continue  # skip header row
                tot_ctr += 1
                l = l.strip().replace(' ', ',')   # covert to CSV
                cols = l.split(',')
                run = int(cols[0])
                subrun = int(cols[1])
                event = int(cols[2])
                energy = float(cols[46])
                for key in pi0_weight_dict:
                    if key[0] == run and key[1] == subrun and key[2] == event and abs(key[3]-energy) < 0.01:
                        wgh_ctr += 1
                        weight = pi0_weight_dict[key]
                        break
                else:
                    if cols[52]=='True':
                        print "WARNING: no weight found for RSE (%s, %s, %s) in file %s! setting to 1..." % (run, subrun, event, josh_list[i])
                        err_ctr += 1
                    else:
                        wgh_ctr += 1  # I think this essentially correct...
                    weight = 1.0
                cols = cols[:54] + cols[56:]  # remove "label" and "filetag" (Josh's legend label and sample label)
                x = ','.join(cols)
                x = x.replace('True', '1')      # convert True/False to 1/0
                x = x.replace('False', '0')     # convert True/False to 1/0
                x = x.replace(',,', ',-9999,')  # convert empty entries to -9999 -- I think this is just where eta_reco is a NaN
                x = x.replace('-3.4028235e+38', '-9999')  # convert what I assume is -inf to -9999 -- I think this is just min/max shower fraction
                if x[-1] != ',': x += ','
                x += "%s\n" % weight
                out.write(x)
    print "INFO: added pi0 weights for %s events" % wgh_ctr
    if ( err_ctr > 0 ):
        print "WARNING: pi0 weights for remaining %s events set to 1" % err_ctr
    if ( wgh_ctr + err_ctr != tot_ctr ):
        print "WARNING: something doesn't math here...!"
    if ( wgh_ctr != len(pi0_weight_dict) ):
        print "WARNING: these values are not consistent!"
        if ( wgh_ctr < len(pi0_weight_dict) ):
            print "  ... but probably fine beacuse events with weights (%s) are less than list of input weights (%s)" % (wgh_ctr, len(pi0_weight_dict))

print "Done!"

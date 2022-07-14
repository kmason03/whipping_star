import os
#Mar24_pi0_weights_1e1p_final_avgscore_ncpi0_run3.txt

events_list  = [ "eventlist_overlay_fset12_run1.txt", "eventlist_intrinsics_fset12_run1.txt", "eventlist_ncpi0_fset12_run1.txt", "eventlist_ccpi0_fset12_run1.txt", "eventlist_ext_fset12_run1.txt",
               "eventlist_overlay_fset12_run2.txt", "eventlist_intrinsics_fset12_run2.txt", 
               "eventlist_overlay_fset12_run3.txt", "eventlist_intrinsics_fset12_run3.txt", "eventlist_ncpi0_fset12_run3.txt", "eventlist_ccpi0_fset12_run3.txt", "eventlist_ext_fset12_run3.txt",
               "eventlist_fullosc_fset12_run1.txt","eventlist_fullosc_fset12_run2.txt","eventlist_fullosc_fset12_run3.txt" ]
weights_list = [ "pi0_weights_1e1p_{}_numu_run1.txt", "pi0_weights_1e1p_{}_nue_run1.txt", "pi0_weights_1e1p_{}_ncpi0_run1.txt", "pi0_weights_1e1p_{}_ccpi0_run1.txt", "empty.txt",
               "pi0_weights_1e1p_{}_numu_run2.txt", "pi0_weights_1e1p_{}_nue_run2.txt", 
               "pi0_weights_1e1p_{}_numu_run3.txt", "pi0_weights_1e1p_{}_nue_run3.txt", "pi0_weights_1e1p_{}_ncpi0_run3.txt", "pi0_weights_1e1p_{}_ccpi0_run3.txt", "empty.txt",
               "pi0_weights_1e1p_{}_fullosc_run1.txt","pi0_weights_1e1p_{}_fullosc_run2.txt","pi0_weights_1e1p_{}_fullosc_run3.txt" ]
out_list   = [ "sel1e1p_{}_bnb_run1_events.txt", "sel1e1p_{}_nue_run1_events.txt", "sel1e1p_{}_ncpi0_run1_events.txt", "sel1e1p_{}_ccpi0_run1_events.txt", "sel1e1p_{}_extbnb_run1_events.txt",
               "sel1e1p_{}_bnb_run2_events.txt", "sel1e1p_{}_nue_run2_events.txt",
               "sel1e1p_{}_bnb_run3_events.txt", "sel1e1p_{}_nue_run3_events.txt", "sel1e1p_{}_ncpi0_run3_events.txt", "sel1e1p_{}_ccpi0_run3_events.txt", "sel1e1p_{}_extbnb_run3_events.txt",
               "sel1e1p_{}_fullosc_run1_events.txt","sel1e1p_{}_fullosc_run2_events.txt","sel1e1p_{}_fullosc_run3_events.txt" ]

sel_list   = [ "FinalSelection" ]
weights_dict = { "FinalSelection": "avgscore"}

top_dir   = "/uboone/app/users/kmason/LaurenVersion/form_input_to_sbnfit"
nick_dir  = os.path.join(top_dir, "rawlists")
weights_dir = os.path.join(top_dir, "fromKatie")

#for s in sel_list:
for s in ["FinalSelection"]:

    print "for selection %s..." % s

    for i in range(len(out_list)):
        
        print "for file %s..." % out_list[i].format(s)

        # Read in the pi0 weights
        pi0_weight_dict = {}
        with open( os.path.join(weights_dir, weights_list[i].format(weights_dict[s])), 'r' ) as k:
            for l in k:
                if l.split(',')[0] == "run": continue
                cols = l.strip().split(',')
                run = int(cols[0])
                subrun = int(cols[1])
                event = int(cols[2])
                # skip true neutrino energy (but good to have it in case we need it) (but we seem to not...?)
                weight = float(cols[4])
                rse = tuple( (run, subrun, event) )
                pi0_weight_dict[rse] = weight

        # Write out the file with the reco vars plus the weight
        tot_ctr = 0
        wgh_ctr = 0
        err_ctr = 0
        this_out = out_list[i].format(s)
        with open( os.path.join(top_dir, this_out), 'w' ) as out:
            nick_fname = events_list[i]
            with open( os.path.join(nick_dir, nick_fname), 'r' ) as n:
                for l in n:
                    if l.startswith('run'): continue  # skip header row
                    tot_ctr += 1
                    l = l.strip().replace(' ', ',')   # covert to CSV
                    cols = l.split(',')
                    run = int(cols[0])
                    subrun = int(cols[1])
                    event = int(cols[2])
                    rse = tuple( (run, subrun, event) )
                    if rse in pi0_weight_dict:
                        wgh_ctr += 1
                        weight = pi0_weight_dict[rse]
                    else:
                        if 'ext' not in out_list[i]:
                            print "WARNING: no weight found for RSE (%s, %s, %s) in file %s! setting to 1..." % (run, subrun, event, nick_list[i])
                        err_ctr += 1
                        weight = 1.0
                    if s == "FinalSelection":
                        cols = cols[:49] + cols[50:]  # remove "label" (Nick's legend label)
                        cols = cols[:52] + cols[53:]  # remove "filetag" (Nick's sample label)
                    else:
                        cols = cols[:49] + cols[50:]  # remove "label" (Nick's legend label)
                        cols = cols[:53] + cols[54:]  # remove "filetag" (Nick's sample label)
                    x = ','.join(cols)
                    x = x.replace('True', '1')      # convert True/False to 1/0
                    x = x.replace('False', '0')     # convert True/False to 1/0
                    x = x.replace(',,', ',-9999,')  # convert empty entries to -9999 -- I think this is just where eta_reco is a NaN
                    if x[-1] != ',': x += ','
                    x += "%s\n" % weight
                    out.write(x)
        print "INFO: added pi0 weights for %s events" % wgh_ctr
        if ( err_ctr > 0 ):
            if 'ext' not in out_list[i]:
                print "WARNING: pi0 weights for remaining %s events set to 1" % err_ctr
        if ( wgh_ctr + err_ctr != tot_ctr ):
            print "WARNING: something doesn't math here...!"
        if ( wgh_ctr != len(pi0_weight_dict) ):
            print "WARNING: these values are not consistent!"

print "Done!"

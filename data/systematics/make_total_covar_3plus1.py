import os,subprocess
import math
import ROOT
from math import sqrt, isnan

def collapseFracCovar(frac_covar, spec, coll_list, debug=False):

    ## Function that calcuates collapsed fractional covariance matrix
    ## Inputs:
    ##   frac_covar    full covariance matrix (TMatrixD)
    ##   spec          spectrum (python list or similar)
    ##   coll_list     information on bins to be collapsed (python list of tuples)
    ##                   note: passed directly to collapseFullCovar and collapseSpec, see there for more info...
    ##   debug         flag for whether to print out various debugging information
    ## Returns:
    ##   coll_covar    collapsed fractional covariance matrix (TMatrixD)

    # Calculate number of bins at input, output
    Nbins_in  = sum([x[0]*x[1] for x in coll_list])
    Nbins_out = sum([x[1] for x in coll_list])
    # If debug, print this out
    if debug:
        print("Number of bins in uncollapsed covariance matrix: {} (should equal {})".format(Nbins_in, frac_covar.GetNrows()))
        print("Number of bins in collaposed covariance matrix: {}".format(Nbins_out))
    # Check dimensions
    err_flag = False
    if Nbins_in != frac_covar.GetNrows():
        print("ERROR: dimensions of frac_covar and coll_list did not match!")
        print(Nbins_in, frac_covar.GetNrows())
        err_flag = True
    if frac_covar.GetNrows() != len(spec):
        print(len(spec),frac_covar.GetNrows())
        print("ERROR: length of spec and coll_list did not match!")
        err_flag = True
    if err_flag:
        return

    # Get the corresponding full covariance matrix
    full_covar = getFullCovar(frac_covar, spec, debug=debug)
    # Collapse the full covariance matrix
    coll_full_covar = collapseFullCovar(full_covar, coll_list, debug=debug)
    # Collapse the spectrum, too
    coll_spec = collapseSpec(spec, coll_list, debug=debug)
    # Get the corresponding collapsed fractional covariance matrix
    coll_frac_covar = getFracCovar(coll_full_covar, coll_spec, debug=debug)

    return coll_frac_covar

def collapseFullCovar(full_covar, coll_list, debug=False):

    ## Function that calcuates collapsed covariance matrix
    ## Inputs:
    ##   full_covar    full covariance matrix (TMatrixD)
    ##   coll_list     information on bins to be collapsed (python list of tuples)
    ##                   ex: if there are three 10-bin channels to collapse followed by a 20-bin channel,
    ##                       then this should be [(3, 10), (1, 20)]
    ##   debug         flag for whether to print out various debugging information
    ## Returns:
    ##   coll_covar    collapsed covariance matrix (TMatrixD)

    # Calculate number of bins at input, output
    Nbins_in  = sum([x[0]*x[1] for x in coll_list])
    Nbins_out = sum([x[1] for x in coll_list])
    # If debug, print this out
    if debug:
        print("Number of bins in full, uncollapsed covariance matrix: {} (should equal {})".format(Nbins_in, full_covar.GetNrows()))
        print("Number of bins in collaposed covariance matrix: {}".format(Nbins_out))
    # Check dimensions
    if Nbins_in != full_covar.GetNrows():
        print("ERROR: dimensions of full_covar and coll_list did not match!")
        return
    # Initialize a matrix of block covariance matrices
    covar_mat = []
    for i in range(len(coll_list)):
        covar_mat.append([])
        for j in range(len(coll_list)):
            covar_mat[i].append( ROOT.TMatrixD(coll_list[i][1], coll_list[j][1]) )
    # And fill it...
    # Loop over matrix of block matrices...
    for i in range(len(coll_list)):
        for j in range(len(coll_list)):
            if debug:
                print("Doing block ({},{})...".format(i,j))
            # Loop over entries within this (i,j)th block...
            for block_i in range(coll_list[i][1]):
                for block_j in range(coll_list[j][1]):
                    if debug:
                        print("  Doing entry ({},{}) within that block...".format(block_i,block_j))
                    # Loop over channels to be summed...
                    for chan_i in range(coll_list[i][0]):
                        for chan_j in range(coll_list[j][0]):
                            full_i = sum([x[0]*x[1] for x in coll_list[:i]]) + chan_i*coll_list[i][1] + block_i
                            full_j = sum([x[0]*x[1] for x in coll_list[:j]]) + chan_j*coll_list[j][1] + block_j
                            if debug:
                                print("    Adding full matrix entry ({},{}) = {}".format(full_i,full_j,full_covar[full_i][full_j]))
                            covar_mat[i][j][block_i][block_j] += full_covar[full_i][full_j]


    # Initialize output collapsed covariance matrix as an empty matrix of the correct dimensions
    coll_covar = ROOT.TMatrixD(Nbins_out, Nbins_out)
    # And fill it...
    # Loop over matrix of block matrices...
    for i in range(len(coll_list)):
        for j in range(len(coll_list)):
            # Loop over entries within this (i,j)th block...
            for block_i in range(coll_list[i][1]):
                for block_j in range(coll_list[j][1]):
                    coll_i = sum([x[1] for x in coll_list[:i]]) + block_i
                    coll_j = sum([x[1] for x in coll_list[:j]]) + block_j
                    coll_covar[coll_i][coll_j] = covar_mat[i][j][block_i][block_j]
                    if debug:
                        s = "Filling collapsed matrix entry ({},{}) with entry ({},{}) from block ({},{})"
                        print(s.format(coll_i,coll_j,block_i,block_j,i,j))
    return coll_covar

def getSpecList(hist_list):

    ## Function that concatenates a list of Th0Ds into a python list
    ## Inputs:
    ##   hist_list    list of spectra to be concatenated (python list of Th0Ds)
    ## Returns:
    ##   spec         concatenated spectrum (python list)

    # Initilize output spectrum
    spec = []

    # Loop over input spectra...
    for hist in hist_list:
        # Loop over bins in this spectrum...
        for i in range(hist.GetNbinsX()):
            spec.append( hist.GetBinContent(i+1) )

    # Return concatenated spectrum
    return spec

def getFullCovar(frac_covar, spec, debug=False):

    ## Function that calcuates full covariance matrix
    ## Inputs:
    ##   frac_covar    fractional covariance matrix (TMatrixD)
    ##   spec          spectrum (python list or similar)
    ##   debug         flag for whether to print out various debugging information
    ## Returns:
    ##   full_covar    full covariance matrix (TMatrixD)

    # If in debug mode, print out information on inputs
    if debug:
        print("spectrum: ", spec)
        print("fractional covariance matrix: ")
        print(frac_covar[0][0])
        #frac_covar.Print()

    # Initialize output full covariance matrix as a copy of input fractional covariance matrix
    full_covar = ROOT.TMatrixD(frac_covar)

    # Compute M_ij = F_ij*N_i*N_j
    for i in range(full_covar.GetNrows()):
        for j in range(full_covar.GetNcols()):
            if isnan(full_covar[i][j]):
                full_covar[i][j] = 0.
            else:
                full_covar[i][j] *= spec[i]*spec[j]

    # If in debug mode, print out information on output
    if debug:
        print("full covariance matrix: ")
        print(full_covar[0][0])
        #full_covar.Print()

    return full_covar

def collapseSpec(full_spec, coll_list, debug=False):

    ## Function that calcuates collapsed spectrum
    ## Inputs:
    ##   full_spec    full spectrum matrix (python list or similar)
    ##   coll_list    information on bins to be collapsed (python list of tuples)
    ##                    ex: if there are three 10-bin channels to collapse followed by a 20-bin channel,
    ##                        then this should be [(3, 10), (1, 20)]
    ##   debug        flag for whether to print out various debugging information
    ## Returns:
    ##   coll_spec    collapsed covariance matrix (TMatrixD)

    # Calculate number of bins at input, output
    Nbins_in  = sum([x[0]*x[1] for x in coll_list])
    Nbins_out = sum([x[1] for x in coll_list])
    # If debug, print this out
    if debug:
        print("Number of bins in uncollapsed spectrum: {} (should equal {})".format(Nbins_in, len(full_spec)))
        print("Number of bins in collaposed spectrum: {}".format(Nbins_out))
    # Check dimensions
    if Nbins_in != len(full_spec):
        print("ERROR: dimensions of full_spec and coll_list did not match!")
        return

    # Initialize output collapsed spectrum as an empty list of correct length
    coll_spec = [ 0. for x in range(Nbins_out) ]
    # And fill it...
    # Loop over blocks of collapsed spectrum...
    for i in range(len(coll_list)):
        if debug:
            print("Doing block {}...".format(i))
        # Loop over bins within that spectrum...
        for block_i in range(coll_list[i][1]):
            coll_i = sum([x[1] for x in coll_list[:i]]) + block_i
            if debug:
                print("  Doing entry {} within that block... global index {}".format(block_i, coll_i))
            # Loop over channels to be summed...
            for chan_i in range(coll_list[i][0]):
                full_i = sum([x[0]*x[1] for x in coll_list[:i]]) + chan_i*coll_list[i][1] + block_i
                coll_spec[coll_i] += full_spec[full_i]
                if debug:
                    print("    Adding full spec entry {} = {}".format(full_i, full_spec[full_i]))

    return coll_spec

def getFracCovar(full_covar, spec, debug=False):

    ## Function that calcuates fractional covariance matrix
    ## Inputs:
    ##   full_covar    full covariance matrix (TMatrixD)
    ##   spec          spectrum (python list or similar)
    ##   debug         flag for whether to print out various debugging information
    ## Returns:
    ##   frac_covar    fractional covariance matrix (TMatrixD)

    # If in debug mode, print out information on inputs
    if debug:
        print("spectrum: ", spec)
        print("full covariance matrix: ")
        print(full_covar[0][0])
        #full_covar.Print()

    # Initialize output fractional covariance matrix as a copy of input full covariance matrix
    frac_covar = ROOT.TMatrixD(full_covar)

    # Compute F_ij = M_ij/(N_i*N_j)
    for i in range(full_covar.GetNrows()):
        for j in range(full_covar.GetNcols()):
            if spec[i]*spec[j] == 0.:
                frac_covar[i][j] = float('nan')
            else:
                frac_covar[i][j] *= 1./(spec[i]*spec[j])

    # If in debug mode, print out information on output
    if debug:
        print("fractional covariance matrix: ")
        print(frac_covar[0][0])
        #frac_covar.Print()

    # Return fractional covariance matrix
    return frac_covar

def main():
    # Declare file names for input files
    # topdir = '/uboone/app/users/yatesla/sbnfit/whipping_star/dllee/forUnblinding/'
    # in_h0_spec_fname = os.path.join(topdir, "h0_v48.SBNspec.root")
    # in_h0_spec_fname = os.path.join(topdir, "h0_v48.SBNspec.root")
    # in_covar_fname   = os.path.join(topdir, "h0_v48.SBNcovar.root")
    fullbins = False
    topdir = ""
    in_h0_spec_fname = ""
    in_covar_fname = ""
    if fullbins:
        topdir='/cluster/tufts/wongjiradlab/jmills09/whipping_star/data/systematics/1xdata/'
        in_h0_spec_fname = os.path.join(topdir,"DL.SBNspec.root")
        in_covar_fname   = os.path.join(topdir,"DL.SBNcovar.root")
    else:
        topdir='/cluster/tufts/wongjiradlab/jmills09/whipping_star/data/systematics/1bin_uncombined_data/'
        in_h0_spec_fname = os.path.join(topdir,"DL_ONEBIN_UNCOMBINED.SBNspec.root")
        in_covar_fname   = os.path.join(topdir,"DL_ONEBIN_UNCOMBINED.SBNcovar.root")

    detsys_fname     = os.path.join(topdir,"fracdetvar_np_19bins.csv")
    # bkg_pred_fname   = os.path.join("/uboone/app/users/yatesla/sbnfit/whipping_star/dllee/bkg/bkg_0.95_prediction.txt")
    # bkg_covar_fname  = os.path.join("/uboone/app/users/yatesla/sbnfit/whipping_star/dllee/bkg/bkg_0.95_cov.txt")

    # Define a few helper variables...
    #   Note: Making some assumptions about the xml configuration
    Nbins_m = 0
    if fullbins:
        Nbins_m = 19
    else:
        Nbins_m = 1
    Nchannels_m = 3

    # offset for each channel (this is messy for combined!)
    # this is the order set in the xml that made the covariance matrix

    offset_1m1p_bnb     = 0
    offset_1m1p_nue     = 1*Nbins_m
    offset_1m1p_ext     = 2*Nbins_m

    # which detector variation each one needs:
    # offset_1m1p_bnb: 1m1p detvar
    # offset_1m1p_nue: none
    # offset_1m1p_extbnb: none


    # Open the input SBNfit files
    in_h0_spec_f = ROOT.TFile.Open(in_h0_spec_fname, "READ")
    in_covar_f   = ROOT.TFile.Open(in_covar_fname, "READ")

    # Declare file names for output files
    if fullbins:
        out_h0_spec_fname = "JOSH_1m1p.SBNspec.root"
        out_covar_fname   = "JOSH_1m1p.SBNcovar.root"
    else:
        out_h0_spec_fname = "JOSH_1m1p_1Bin.SBNspec.root"
        out_covar_fname   = "JOSH_1m1p_1Bin.SBNcovar.root"

    # Create the output files
    out_h0_spec_f = ROOT.TFile(out_h0_spec_fname, "RECREATE")
    out_covar_f   = ROOT.TFile(out_covar_fname, "RECREATE")

    # Initialize the output spectra, covariance matrix as copies of the inputs
    out_h0_spec_dict = {}
    for k in [ key.GetName() for key in in_h0_spec_f.GetListOfKeys() ]:
        out_h0_spec_dict[k] = ROOT.TH1D( in_h0_spec_f.Get(k) )
    out_covar = ROOT.TMatrixD( in_covar_f.Get("frac_covariance") )
    # Zero out any nans, so that things we add actually get added
    for i in range(out_covar.GetNrows()):
        for j in range(out_covar.GetNcols()):
            if math.isnan(out_covar[i][j]):
                out_covar[i][j] = 0.

    # Add detector systematics for 1e1p nue, 1e1p lee, and 1m1p bnb...
    # Read in fractional detector systematic covariance matrix from csv file
    #   Note: Covariance matrix in this file has dimension 20+12
    #           20 1m1p bins from 200 to 1200 MeV in  50 MeV bins, then
    #           22 1e1p  bins from   200 to 2400 MeV in 100 MeV bins
    print("Adding detector systematics from {}".format(detsys_fname))
    detsys_covar = []
    with open(detsys_fname, 'r') as f:
        for l in f:
            detsys_covar.append( [ float(x) for x in l.strip().split(',') ] )

    # Break this out into block matrices with correct dimensions - the one from Ran is backwards
    #   Note: We only use the 19 1m1p bins from 250 to 1200 MeV, so add 1 to the offset
    # these are the offsets of where blocks are in the detvar matrix, not to be confused with the total

    # 1bin hack
    detsys_covar_mm = [[0]]
    if not fullbins:
        for i in range(19):
            for j in range(19):
                detsys_covar_mm[0][0] += detsys_covar[i][j]
    else:
        detsys_covar_mm = detsys_covar




    # for i in range(Nbins_m):
    #     detsys_covar_mm.append( [] )
    #     for j in range(Nbins_m):
    #         detsys_covar_mm[i].append( detsys_covar[in_offset_m+i][in_offset_m+j] )





    # Add detsys_covar_mm to the 1m1p bnb,1m1p ccpi0, 1m1p ncpi0
    for i in range(Nbins_m):
        for j in range(Nbins_m):
            # main block
            out_covar[offset_1m1p_bnb+i][offset_1m1p_bnb+j] += detsys_covar_mm[i][j]

    # ... done adding detector systematics for 1e1p nue, 1e1p fullosc, 1m1p bnb, 1m1p ccpi0, 1m1p ncpi0

    # Add *fractional* mc stat errors for all subchannels
    #   and remove them from spec errors, so SBNfit doesn't double-count


    # For Stat Only
    # for i in range(out_covar.GetNrows()):
    #     for j in range(out_covar.GetNcols()):
    #         out_covar[i][j] = 0.


    print("Adding mc stat errors from {}".format(in_h0_spec_fname))
    # first the 1e1p (doing it concurrently)
    # next is the 1m1p
    for i in range(Nbins_m):
        # bnb
        if in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinContent(i+1) > 0.:
            print(in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinContent(i+1), in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinError(i+1), in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinError(i+1)/ in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinContent(i+1))
            out_covar[i+offset_1m1p_bnb][i+offset_1m1p_bnb] += (in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinError(i+1) / in_h0_spec_f.Get("nu_uBooNE_1m1p_bnb").GetBinContent(i+1))**2
        out_h0_spec_dict["nu_uBooNE_1m1p_bnb"].SetBinError(i+1, 0.)
        # nue
        if in_h0_spec_f.Get("nu_uBooNE_1m1p_nue").GetBinContent(i+1) > 0.:
            out_covar[i+offset_1m1p_nue][i+offset_1m1p_nue] += (in_h0_spec_f.Get("nu_uBooNE_1m1p_nue").GetBinError(i+1) / in_h0_spec_f.Get("nu_uBooNE_1m1p_nue").GetBinContent(i+1))**2
        out_h0_spec_dict["nu_uBooNE_1m1p_nue"].SetBinError(i+1, 0.)
        # ext
        if in_h0_spec_f.Get("nu_uBooNE_1m1p_ext").GetBinContent(i+1) > 0.:
            out_covar[i+offset_1m1p_ext][i+offset_1m1p_ext] += (in_h0_spec_f.Get("nu_uBooNE_1m1p_ext").GetBinError(i+1) / in_h0_spec_f.Get("nu_uBooNE_1m1p_ext").GetBinContent(i+1))**2
        out_h0_spec_dict["nu_uBooNE_1m1p_ext"].SetBinError(i+1, 0.)

    # ... done adding mc stat errors for 1e1p nue, 1e1p lee, and 1m1p

    # Write everything out
    for k in [ key.GetName() for key in in_h0_spec_f.GetListOfKeys() ]:
        out_h0_spec_f.WriteTObject(out_h0_spec_dict[k])
    out_covar_f.WriteTObject(out_covar, "frac_covariance")
    # next get the collapsed covariance matrix
    joint_keys    = [ key.GetName() for key in in_h0_spec_f.GetListOfKeys() ]
    joint_spec = getSpecList( [ out_h0_spec_dict[k] for k in joint_keys ] )
    # make sure to change the binning here for your subchannels
    out_covar_collapse = 0
    if fullbins:
        out_covar_collapse = collapseFracCovar(out_covar, joint_spec, [(3,19)], True)
    else:
        out_covar_collapse = collapseFracCovar(out_covar, joint_spec, [(3,1)], True)

    out_covar_f.WriteTObject(out_covar_collapse, "frac_covariance_collapsed")
    # Close the outputs
    out_h0_spec_f.Close()
    out_covar_f.Close()
    # Close the inputs
    in_h0_spec_f.Close()
    in_covar_f.Close()

    print("Done!")

if __name__ == "__main__":
    main()

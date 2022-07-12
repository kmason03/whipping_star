#!/bin/bash

source /usr/local/root/root-6.22.06/bin/thisroot.sh
makewhip()
{
        pwd
	cd /cluster/tufts/wongjiradlab/jmills09/whipping_star/build/
	cmake ../
	make
	#cd DecayNuMuDis/
	cd NuMuDisappearance/
}

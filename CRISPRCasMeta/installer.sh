#!/bin/bash

#------------------------
function launchInstall {
# $1 $packageManagmentInstall
# $2 package name
# $3 $LOGFILE
    echo "Installation of $2" >> $3
    $1 $2 >> $3
}
function checkProcess {
# $1 Step
# $2 exit error level
	if [ $? -ne 0 ];
    then
        echo "$1 : FAIL"
        exit $2
	else
		echo "$1 : OK"
    fi
}
#------------------------

sudo chmod 755 .
sudo chmod 755 *

CURDIR=`pwd` 
LOGFILE=$CURDIR/install.log
if [ -e $LOGFILE ]; then rm $LOGFILE ; fi

ostype=`echo $OSTYPE`
echo "$ostype" >> $LOGFILE

if [ ! "$ostype" = "linux-gnu" ]; then
    echo 'Sorry, install process only for OSTYPE linux-gnu (based on apt and apt-get)'
else
    ##Python library and make command
    echo ""
    echo "step (1/4) Installation of python regex library and command make..."
    echo ""
	apt-get update
    #sudo apt -y install python-pip >> $LOGFILE #try install python or update it
    #sudo pip install regex >> $LOGFILE
    type pip || sudo apt -y install python-pip >> $LOGFILE
    checkProcess "Python pip command test/installation " 3
    sudo pip install regex >> $LOGFILE
    checkProcess "installation of regex library (Python 2.7) " 4
    #sudo apt -y install make >> $LOGFILE
    type make || sudo apt -y install make >> $LOGFILE
    checkProcess "command make test/installation " 3

    ##CRT(Minced)
    echo ""
    echo "step (2/4) Installation of CRT(Minced)..."
    echo ""
    #sudo apt -y install default-jdk >> $LOGFILE
	if [ ! -d minced-master/ ]; then
		type javac || sudo apt -y install default-jdk >> $LOGFILE
		checkProcess "Java default jdk test/installation " 3
		#wget https://github.com/ctSkennerton/minced/archive/master.zip
		wget https://github.com/ctSkennerton/minced/archive/0.3.0.zip
		#unzip master.zip >> $LOGFILE
		unzip 0.3.0.zip >> $LOGFILE
		#cd minced-master/
		cd minced-0.3.0/
		make >> $LOGFILE
		checkProcess "CRT(MinCED) installation " 3
		cd ..
	fi

    ##CRISPRCasFinder
    echo ""
    echo "step (3/4) Installation of MetaFinder (CRISPRCasFinder metagenomic version)..."
    echo ""
    packageManagmentInstall='sudo apt-get -y install '
    #distribution='Linux_x86_64'
    launchInstall "$packageManagmentInstall" "cpanminus" "$LOGFILE"
    sudo cpanm -f Bio::Tools::Run::Alignment::Muscle >> $LOGFILE
    checkProcess "Installation of Bio::Tools::Run::Alignment::Muscle library (Perl) " 4
    #sudo apt-get -y install npm >> $LOGFILE
    type npm || sudo apt-get -y install npm >> $LOGFILE
    sudo npm install parse-json >> $LOGFILE
    sudo cpanm --force JSON::Parse >> $LOGFILE
    checkProcess "Installation of JSON::Parse (Perl) " 4
    cd MetaFinder/
	echo "You will have to press ENTER key if it seems to take a long time during next steps (more than 5 minutes for a step) !"
    bash installer_UBUNTU.sh >> $LOGFILE
    checkProcess "MetaFinder, installation " 3
    source ~/.profile
    perl CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def General -cas -i install_test/sequence.fasta -out Results_test_install –keep
    cd ..

    #Test CRISPRCasMeta
    echo ""
    echo "step (4/4) CRISPRCasMeta installation test..."
    echo ""
    python CRISPRCasMeta.py -i $CURDIR/DataTest/B2F01_1000.fasta.gz -o Test_Output
    echo "CRISPRCasMeta, difference with expected results : "
    diff $CURDIR/DataTest/Test_Output/result.json $CURDIR/DataTest/result.json
    #Test CRISPRCasFinder
    cd MetaFinder/
    echo "CRISPRCasFinder, difference with expected results : "
    diff Results_test_install/TSV/Cas_REPORT.tsv install_test/Cas_REPORT.tsv
    diff Results_test_install/TSV/Crisprs_REPORT.tsv install_test/Crisprs_REPORT.tsv

fi 
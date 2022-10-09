#!/bin/bash
## Adjusted version of the script by Ben Bolker,
## https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/misc/autowinbuilder.sh
## 
## For the main package maintainer, use
##   devtools::check_win_release()
##   devtools::check_win_devel()
## instead of this script!
## 
## This script temorarily changes the maintainer email of current branch
## to finn.lindgren@gmail.com,
## builds the tarball, and uploads to win-builder
## run from the main excursions folder, as
##   misc/autowinbuilder.sh
## 
## arg 1: which platform to test (both, release, or devel) ?

whichrel=${1:-both}
echo $whichrel
MY_EMAIL=finn.lindgren@gmail.com
MAINTAINER_EMAIL=davidbolin@gmail.com
version=`grep 'Version' DESCRIPTION | sed -e 's/Version: //'`
echo "excursions version $version"
tarball="excursions_${version}.tar.gz"
## substitute e-mail in DESCRIPTION file and build tarball
sed -i -e "s/$MAINTAINER_EMAIL/$MY_EMAIL/" DESCRIPTION
(cd .. && R CMD build excursions) || exit 1
tar zxvfO ../$tarball excursions/DESCRIPTION | grep Maintainer
echo "tarball in ..: $tarball"
## https://serverfault.com/questions/279176/ftp-uploading-in-bash-script
HOST=win-builder.r-project.org
USER=ftp
PASS=$MY_EMAIL
if [ $whichrel == "both" ] || [ $whichrel == "devel" ]; then
echo "uploading to win-builder/devel"
(cd .. && ftp -inv $HOST << EOT

user $USER $PASS
binary
cd R-devel
put $tarball
bye 
EOT
)
fi

if [ $whichrel == "both" ] || [ $whichrel == "release" ]; then
echo "uploading to win-builder/release"
## upload to R-release
(cd .. && ftp -inv $HOST << EOT

user $USER $PASS
binary
cd R-release
put $tarball
bye
EOT
)
fi

## revert changes to DESCRIPTION file
git checkout -- DESCRIPTION

## check status
Rscript -e "foghorn::winbuilder_queue('excursions')"

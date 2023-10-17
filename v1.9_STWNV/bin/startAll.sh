#!/bin/bash
set -e

cd ../v1.9_s1/bin
source start.sh

cd ../../v1.9_s2/bin
source start.sh

cd ../../v1.9_s3/bin
source start.sh

cd ../../v1.9_s4/bin
source start.sh

cd ../../v1.9_s5/bin
source start.sh

echo "========================================="
echo "All Cases from s1 to s5 started"
echo "========================================="

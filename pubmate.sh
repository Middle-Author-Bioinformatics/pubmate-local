#!/bin/bash
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

exec > >(tee -i /home/ark/MAB/pubmate/pubmate_looper.log)
exec 2>&1

## Debugging information
#echo "Script started at $(date)"
#echo "Current directory: $(pwd)"
#echo "Environment variables: $(env)"

eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate base  # Activate the base environment where `boto3` is installed

# Debug PATH
#echo "Updated PATH: $PATH"

# Debug Python environment
#echo "Python version being used:"
#which python3
#python3 --version
#
#echo "Python modules installed:"
#python3 -m pip list

KEY=$1
ID=$KEY
DIR=/home/ark/MAB/pubmate/${ID}
OUT=/home/ark/MAB/pubmate/completed/${ID}-results

mkdir ${OUT}

name=$(grep 'Name' ${DIR}/form-data.txt | cut -d ' ' -f2)
email=$(grep 'Email' ${DIR}/form-data.txt | cut -d ' ' -f2)
task=$(grep 'Tab' ${DIR}/form-data.txt | cut -d ' ' -f2)

# Verify email
result=$(python3 /home/ark/MAB/bin/pubmate-local/check_email.py --email ${email})
echo $result

# Set PATH to include Conda and script locations
export PATH="/home/ark/miniconda3/bin:/usr/local/bin:/usr/bin:/bin:/home/ark/MAB/bin/pubmate-local:$PATH"
eval "$(/home/ark/miniconda3/bin/conda shell.bash hook)"
conda activate pubmate

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate Conda environment."
    exit 1
fi
sleep 5


if [[ ${task} == "topic" ]]; then
    keywords=$(grep 'Keywords' ${DIR}/form-data.txt | cut -d ' ' -f2)
    question=$(grep 'Question' ${DIR}/form-data.txt | cut -d ' ' -f2-)
    echo /home/ark/MAB/bin/pubmate-local/GetTheGist.py --keywords "${keywords}" --question "${question}" --output1 ${OUT}/abstracts.pdf ----output2 ${OUT}/gpt_says.pdf
    /home/ark/MAB/bin/pubmate-local/GetTheGist.py --keywords "${keywords}" --question "${question}" --output1 ${OUT}/abstracts.pdf --output2 ${OUT}/gpt_says.pdf
elif [[ ${task} == 'author' ]]; then
    first=$(grep 'First' ${DIR}/form-data.txt | cut -d ' ' -f2)
    middle=$(grep 'Middle' ${DIR}/form-data.txt | cut -d ' ' -f2)
    last=$(grep 'Last' ${DIR}/form-data.txt | cut -d ' ' -f2)
    topic=$(grep 'Topic' ${DIR}/form-data.txt | cut -d ' ' -f2-)
    echo /home/ark/MAB/bin/pubmate-local/pubcard.py --first "${first}" --middle "${middle}" --last "${last}" --topic "${topic}"--outputImage ${OUT}/pubcard.png --outputStats ${OUT}/pubcard.pdf --output ${OUT}/papers.pdf
    /home/ark/MAB/bin/pubmate-local/pubcard.py --first "${first}" --middle "${middle}" --last "${last}" --topic "${topic}"--outputImage ${OUT}/pubcard.png --outputStats ${OUT}/pubcard.pdf --output ${OUT}/papers.pdf
else
    echo "--- Invalid task: ${task}"
fi

# **************************************************************************************************
# **************************************************************************************************
# **************************************************************************************************
if [ $? -ne 0 ]; then
    echo "Error: PubMate failed."
    conda deactivate
    exit 1
fi
conda deactivate
sleep 5

# Archive results
mv /home/ark/MAB/pubmate/completed/${ID}-results ./${ID}-results
tar -cf ${ID}-results.tar ${ID}-results && gzip ${ID}-results.tar

# Upload results to S3 and generate presigned URL
results_tar="${ID}-results.tar.gz"
s3_key="${ID}-results.tar.gz"
python3 /home/ark/MAB/bin/pubmate-local/push.py --bucket binfo-dump --output_key ${s3_key} --source ${results_tar}
url=$(python3 /home/ark/MAB/bin/pubmate-local/gen_presign_url.py --bucket binfo-dump --key ${s3_key} --expiration 86400)

mv ${ID}-results.tar.gz /home/ark/MAB/pubmate/completed/${ID}-results.tar.gz

# Send email
python3 /home/ark/MAB/bin/pubmate-local/send_email.py \
    --sender ark@midauthorbio.com \
    --recipient ${email} \
    --subject "Your PubMate Results!" \
    --body "Hi ${name},

    Your PubMate results are available for download using the link below. The link will expire in 24 hours.

    ${url}

    Thanks!
    Ark"

if [ $? -ne 0 ]; then
    echo "Error: send_email.py failed."
    conda deactivate
    exit 1
fi

sleep 5

#sudo rm -rf ${DIR}

conda deactivate
echo "PubMate completed successfully."




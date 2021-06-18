echo "running unit tests:"
python test/unit/samplot_test.py
echo "finished unit tests"
echo "running functional tests for \`plot\`:"
bash test/func/samplot_test.sh
printf "\n\nfinished functional tests for \`plot\`:\n"
printf "running functional tests for \`vcf\`:\n"
bash test/func/samplot_vcf_test.sh
echo "finished functional tests for \`vcf\`:"

###########
metab5e_8="cojoinput/metab"
genome="genome1-22/"
output_dir="cojoresult/metab/"

# 创建一个函数，用于执行单个文件的分析
perform_analysis() {
    file="$1"
    FILE=$(echo "$file" | awk -F'.' '{print $2}')
    echo "$FILE"
    mkdir -p ${output_dir}/${FILE}_clumpped;

    for m in {1..22}; do
        gcta-1.94.1 --bfile ${genome}/genochr${m} --chr ${m} --maf 0.01 --cojo-p 5e-8 --cojo-wind 5000 --cojo-collinear 0.9 --cojo-slct --thread-num 16 --cojo-file "$file" --out ${output_dir}/${FILE}_clumpped/${FILE}_chr"${m}"
        if [ $? -ne 0 ]; then
            echo "Error in GCTA analysis for ${FILE}_chr${m}."
            exit 1
        fi
    done;
}

# 并行执行脚本
for file in ${metab5e_8}/*; do
    perform_analysis "$file" &

    # 控制并发任务数量
    if (( $(jobs | wc -l) >= 22 )); then
        wait  # 等待前面的并发任务完成
    fi
done

# 等待剩余的并发任务完成
wait


###########

metab5e_8="cojoinput/ratios"
genome="genome1-22/"
output_dir="cojoresult/ratios/"

# 创建一个函数，用于执行单个文件的分析
perform_analysis() {
    file="$1"
    FILE=$(echo "$file" | awk -F'.' '{print $2"."$3}')
    echo "$FILE"
    mkdir -p ${output_dir}/${FILE}_clumpped;

    for m in {1..22}; do
        gcta-1.94.1 --bfile ${genome}/genochr${m} --chr ${m} --maf 0.01 --cojo-p 5e-8 --cojo-wind 5000 --cojo-collinear 0.9 --cojo-slct --thread-num 16 --cojo-file "$file" --out ${output_dir}/${FILE}_clumpped/${FILE}_chr"${m}"
        if [ $? -ne 0 ]; then
            echo "Error in GCTA analysis for ${FILE}_chr${m}."
            exit 1
        fi
    done;
}

# 并行执行脚本
for file in ${metab5e_8}/*; do
    perform_analysis "$file" &

    # 控制并发任务数量
    if (( $(jobs | wc -l) >= 30 )); then
        wait  # 等待前面的并发任务完成
    fi
done

# 等待剩余的并发任务完成
wait




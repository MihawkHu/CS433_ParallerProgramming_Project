./bin/hdfs dfs -rm /input/init.txt
./bin/hdfs dfs -put nbody/init.txt /input/
./bin/hdfs dfs -rm -r out
./bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file nbody/mapper.py -mapper "python nbody/mapper.py" -file nbody/reducer.py -reducer "python nbody/reducer.py" -input /input/init.txt  -output out
rm -r out
./bin/hdfs dfs -get out out
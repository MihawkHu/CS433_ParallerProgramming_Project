./bin/hdfs dfs -rm -r out
./bin/hadoop jar share/hadoop/tools/lib/hadoop-streaming-2.7.3.jar -file nbody/mapper.py -mapper "python InvertedIndex/mapper.py" -file nbody/reducer.py -reducer "python InvertedIndex/reducer.py" -input /input  -output out
rm -r out
./bin/hdfs dfs -get out out

# Coreset for Individual Fairness

Compile the code using:
```
make
```

## Run Individual Fairness Coreset Algorithm
```
./mainIF_Coreset input_coreset_file_name no_of_centers full_data_file_name fair_radius_file_name_MV
```

### Parameters:
- `input_coreset_file_name`: Input dataset, where each line is a data point -- this is the coreset file or full data file
- `no_of_centers`: Number of centers (k)
- `full_data_file_name`: File containing the full dataset 
- `fair_radius_file_name_MV`: file containing fair radius for each data point on each line. Each line contains two space separated values. First value is an integer which is the ID of the data point, and second value is a float, which is the fair radius according to [Mahabadi Vakilian Algorithm](https://arxiv.org/pdf/2002.06742.pdf).


## Creating coreset from data:
```
python3 createCoreset.py datasetName input_file no_of_centers real_or_semi-synthetic
```
### Parameters:
- `datasetName`: dataset name *(make sure that a folder named 'datasetNameCoreset' exists inside 'datasets' folder)*
- `input_file`: name of the input dataset file
- `no_of_centers`: No. of centers (k)
- `real_or_semi-synthetic`: 0 for real dataset and 1 for semi-synthetic


## Run [Mahabadi Vakilian Algorithm](https://arxiv.org/pdf/2002.06742.pdf):

 ```
./mainIF_MV input_file_name no_of_centers
 ```

## Creating Semi-Synthetic Data:
```
python3 createSemiSyntheticDataset.py input_file_name output_file_name no_of_uar_points no_of_final_generated_points
```
### Parameters:
- `input_file_name`: Input dataset, where each line is a data point
- `output_file_name`: name of the file that will contain the final generated set of points
- `no_of_uar_points`: number of points to be selected uniformly at random (these points will be duplicated using a power law distribution on these points)
- `no_of_final_generated_points`: number of points to be generated at the end using the power law distribution on the uniformly selected points.   

## Creating coreset from Semi-Synthetic Data:
```
python3 createCoreset.py datasetName input_file no_of_centers real_or_semi-synthetic
```
### Parameters:
- `datasetName`: dataset name *(make sure that a folder named 'datasetNameCoreset' exists inside 'datasets' folder)*
- `input_file`: name of the input dataset file
- `no_of_centers`: No. of centers (k)
- `real_or_semi-synthetic`: 1 for semi-synthetic and 0 for real dataset

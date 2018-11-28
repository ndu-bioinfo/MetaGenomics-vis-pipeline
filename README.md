## Functions for creating diversity plots using plotly

#### alpha or beta diversity bar plot

The data frame contains alphda diversity values should be in the following format;  diversity values calculated by different methodolgies can be stored in separated columns and called in the plot funcitons. 

| index	| sample_code	| heip_e	| season	| year	| region|
| --- | --- | --- | --- | --- | --- |
| 201501_086.7_110.0_86	| ND0184	| 0.113271	| Winter	| 2015	| Off|
| 201402_090.0_090.0_10	| ND0012	| 0.103633	| Winter	| 2014	| Off|
| 201607_093.3_110.0_17	| ND0398	| 0.065743	| Summer	| 2016	| Off|
| 201604_093.3_026.7_20	| ND0355	| 0.179908	| Spring	| 2016	| SCB|
| 201407_086.7_045.0_53	| ND0094	| 0.133548	| Summer	| 2014	| Up|
| ...|

beta diversity data are stored in dictionary, with keys as the name of methods used and items as the similarity/distance matrics associated.

|     | 201501_086.7_110.0_86	| 201402_090.0_090.0_10 |	201607_093.3_110.0_17 |	201604_093.3_026.7_20 |	201407_086.7_045.0_53|
| --- | --- | --- | --- | --- | --- |
| 201501_086.7_110.0_86	| 0.000000	| 0.518750	| 0.624583	| 0.781583	| 0.651833| 
| 201402_090.0_090.0_10	| 0.518750	| 0.000000	| 0.480417	| 0.641250	| 0.536583| 
| 201607_093.3_110.0_17	| 0.624583	| 0.480417	| 0.000000	| 0.796750	| 0.761083| 
| 201604_093.3_026.7_20	| 0.781583	| 0.641250	| 0.796750	| 0.000000	| 0.566000| 
| 201407_086.7_045.0_53	| 0.651833	| 0.536583	| 0.761083	| 0.566000	| 0.000000| 



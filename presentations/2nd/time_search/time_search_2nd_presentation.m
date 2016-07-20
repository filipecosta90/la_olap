
red = [ 255 0 0 ];
green = [ 0 255 0 ];
blue = [ 0 0 255 ];
%yellow = [ 255 255 0 ];

yellow = [ 0 128 128 ];

magenta = [ 255 0 255 ];
brown = [ 165 42 42 ];
orange = [ 255 165 0 ];

color0 = red/255;
color1 = green/255;
color2 = blue/255;
color6 = yellow/255;
color7 = magenta/255;
color8 = brown/255;

color9 = orange/255;



FigHandle = figure;
  set(FigHandle, 'Position', [0, 0, 760, 560]);
cd ('../../../src');

dataset_1 = csvread('timing/timings_vec_1.csv');
dataset_2 = csvread('timing/timings_vec_2.csv');
dataset_4 = csvread('timing/timings_vec_4.csv');
dataset_8 = csvread('timing/timings_vec_8.csv');
dataset_16 = csvread('timing/timings_vec_16.csv');
dataset_32 = csvread('timing/timings_vec_32.csv');


time_1 = dataset_1 ( :, 2); 
percntiles_1 = prctile(time_1,[5 95]); %5th and 95th percentile
outlierIndex_1 = time_1 < percntiles_1(1) | time_1 > percntiles_1(2);
%remove outlier values
time_1(outlierIndex_1) = [];
best_time_1 = min( time_1 );

time_2 = dataset_2 ( :, 2); 
percntiles_2 = prctile(time_2,[5 95]); %5th and 95th percentile
outlierIndex_2 = time_2 < percntiles_2(1) | time_2 > percntiles_2(2);
%remove outlier values
time_2(outlierIndex_2) = [];
best_time_2 = min( time_2 );

time_4 = dataset_4 ( :, 2); 
percntiles_4 = prctile(time_4,[5 95]); %5th and 95th percentile
outlierIndex_4 = time_4 < percntiles_4(1) | time_4 > percntiles_4(2);
%remove outlier values
time_4(outlierIndex_4) = [];
best_time_4 = min( time_4 );

time_8 = dataset_8 ( :, 2); 
percntiles_8 = prctile(time_8,[5 95]); %5th and 95th percentile
outlierIndex_8 = time_8 < percntiles_8(1) | time_8 > percntiles_8(2);
%remove outlier values
time_8(outlierIndex_8) = [];
best_time_8 = min( time_8 );

time_16 = dataset_16 ( :, 2); 
percntiles_16 = prctile(time_16,[5 95]); %5th and 95th percentile
outlierIndex_16 = time_16 < percntiles_16(1) | time_16 > percntiles_16(2);
%remove outlier values
time_16(outlierIndex_16) = [];
best_time_16 = min( time_16 );

time_32 = dataset_32 ( :, 2); 
percntiles_32 = prctile(time_32,[5 95]); %5th and 95th percentile
outlierIndex_32 = time_32 < percntiles_32(1) | time_32 > percntiles_32(2);
%remove outlier values
time_32(outlierIndex_32) = [];
best_time_32 = min( time_32 );

seq_pgres_1 = csvread('timing_sql_seq/seq_pgres_1.csv') / 1000;
seq_pgres_2 = csvread('timing_sql_seq/seq_pgres_2.csv') / 1000;
seq_pgres_4 = csvread('timing_sql_seq/seq_pgres_4.csv') / 1000;
seq_pgres_8 = csvread('timing_sql_seq/seq_pgres_8.csv') / 1000;
seq_pgres_16 = csvread('timing_sql_seq/seq_pgres_16.csv') / 1000;
seq_pgres_32 = csvread('timing_sql_seq/seq_pgres_32.csv') / 1000;


percntiles_pgres_1 = prctile(seq_pgres_1,[5 95]); %5th and 95th percentile
outlierIndex_pgres_1 = seq_pgres_1 < percntiles_pgres_1(1) | seq_pgres_1 > percntiles_pgres_1(2);
%remove outlier values
seq_pgres_1(outlierIndex_pgres_1) = [];
best_time_pgres_1 = min( seq_pgres_1 );

percntiles_pgres_2 = prctile(seq_pgres_2,[5 95]); %5th and 95th percentile
outlierIndex_pgres_2 = seq_pgres_2 < percntiles_pgres_2(1) | seq_pgres_2 > percntiles_pgres_2(2);
%remove outlier values
seq_pgres_2(outlierIndex_pgres_2) = [];
best_time_pgres_2 = min( seq_pgres_2 );

percntiles_pgres_4 = prctile(seq_pgres_4,[5 95]); %5th and 95th percentile
outlierIndex_pgres_4 = seq_pgres_4 < percntiles_pgres_4(1) | seq_pgres_4 > percntiles_pgres_4(2);
%remove outlier values
seq_pgres_4(outlierIndex_pgres_4) = [];
best_time_pgres_4 = min( seq_pgres_4 );

percntiles_pgres_8 = prctile(seq_pgres_8,[5 95]); %5th and 95th percentile
outlierIndex_pgres_8 = seq_pgres_8 < percntiles_pgres_8(1) | seq_pgres_8 > percntiles_pgres_8(2);
%remove outlier values
seq_pgres_8(outlierIndex_pgres_8) = [];
best_time_pgres_8 = min( seq_pgres_8 );

percntiles_pgres_16 = prctile(seq_pgres_16,[5 95]); %5th and 95th percentile
outlierIndex_pgres_16 = seq_pgres_16 < percntiles_pgres_16(1) | seq_pgres_16 > percntiles_pgres_16(2);
%remove outlier values
seq_pgres_16(outlierIndex_pgres_16) = [];
best_time_pgres_16 = min( seq_pgres_16 );

percntiles_pgres_32 = prctile(seq_pgres_32,[5 95]); %5th and 95th percentile
outlierIndex_pgres_32 = seq_pgres_32 < percntiles_pgres_32(1) | seq_pgres_32 > percntiles_pgres_32(2);
%remove outlier values
seq_pgres_32(outlierIndex_pgres_32) = [];
best_time_pgres_32 = min( seq_pgres_32 );

par_pgres_1 = csvread('timing_sql_par/par_pgres_1.csv') / 1000;
par_pgres_2 = csvread('timing_sql_par/par_pgres_2.csv') / 1000;
par_pgres_4 = csvread('timing_sql_par/par_pgres_4.csv') / 1000;
par_pgres_8 = csvread('timing_sql_par/par_pgres_8.csv') / 1000;
par_pgres_16 = csvread('timing_sql_par/par_pgres_16.csv') / 1000;
par_pgres_32 = csvread('timing_sql_par/par_pgres_32.csv') / 1000;

percntiles_pgres_par_1 = prctile(par_pgres_1,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_1 = par_pgres_1 < percntiles_pgres_par_1(1) | par_pgres_1 > percntiles_pgres_par_1(2);
%remove outlier values
par_pgres_1(outlierIndex_pgres_par_1) = [];
best_time_pgres_par_1 = min( par_pgres_1 );

percntiles_pgres_par_2 = prctile(par_pgres_2,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_2 = par_pgres_2 < percntiles_pgres_par_2(1) | par_pgres_2 > percntiles_pgres_par_2(2);
%remove outlier values
par_pgres_2(outlierIndex_pgres_par_2) = [];
best_time_pgres_par_2 = min( par_pgres_2 );

percntiles_pgres_par_4 = prctile(par_pgres_4,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_4 = par_pgres_4 < percntiles_pgres_par_4(1) | par_pgres_4 > percntiles_pgres_par_4(2);
%remove outlier values
par_pgres_4(outlierIndex_pgres_par_4) = [];
best_time_pgres_par_4 = min( par_pgres_4 );

percntiles_pgres_par_8 = prctile(par_pgres_8,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_8 = par_pgres_8 < percntiles_pgres_par_8(1) | par_pgres_8 > percntiles_pgres_par_8(2);
%remove outlier values
par_pgres_8(outlierIndex_pgres_par_8) = [];
best_time_pgres_par_8 = min( par_pgres_8 );

percntiles_pgres_par_16 = prctile(par_pgres_16,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_16 = par_pgres_16 < percntiles_pgres_par_16(1) | par_pgres_16 > percntiles_pgres_par_16(2);
%remove outlier values
par_pgres_16(outlierIndex_pgres_par_16) = [];
best_time_pgres_par_16 = min( par_pgres_16 );

percntiles_pgres_par_32 = prctile(par_pgres_32,[5 95]); %5th and 95th percentile
outlierIndex_pgres_par_32 = par_pgres_32 < percntiles_pgres_par_32(1) | par_pgres_32 > percntiles_pgres_par_32(2);
%remove outlier values
par_pgres_32(outlierIndex_pgres_par_32) = [];
best_time_pgres_par_32 = min( par_pgres_32 );

novec_1 =  csvread('timing/timings_no_vec_1.dat') ;
novec_1 =  novec_1( :, 2);
novec_2 = csvread('timing/timings_no_vec_2.dat');
novec_2 =  novec_2( :, 2);
novec_4 = csvread('timing/timings_no_vec_4.dat');
novec_4 =  novec_4( :, 2);
novec_8 = csvread('timing/timings_no_vec_8.dat');
novec_8 =  novec_8( :, 2);
novec_16 = csvread('timing/timings_no_vec_16.dat');
novec_16 =  novec_16( :, 2);
novec_32 = csvread('timing/timings_no_vec_32.dat');
novec_32 =  novec_32( :, 2);


percntiles_novec_1 = prctile(novec_1,[5 95]); %5th and 95th percentile
outlierIndex_novec_1 = novec_1 < percntiles_novec_1(1) | novec_1 > percntiles_novec_1(2);
%remove outlier values
novec_1(outlierIndex_novec_1) = [];
best_time_novec_1 = max( novec_1 );

percntiles_novec_2 = prctile(novec_2,[5 95]); %5th and 95th percentile
outlierIndex_novec_2 = novec_2 < percntiles_novec_2(1) | novec_2 > percntiles_novec_2(2);
%remove outlier values
novec_2(outlierIndex_novec_2) = [];
best_time_novec_2 = max( novec_2 );

percntiles_novec_4 = prctile(novec_4,[5 95]); %5th and 95th percentile
outlierIndex_novec_4 = novec_4 < percntiles_novec_4(1) | novec_4 > percntiles_novec_4(2);
%remove outlier values
novec_4(outlierIndex_novec_4) = [];
best_time_novec_4 = max( novec_4 );

percntiles_novec_8 = prctile(novec_8,[5 95]); %5th and 95th percentile
outlierIndex_novec_8 = novec_8 < percntiles_novec_8(1) | novec_8 > percntiles_novec_8(2);
%remove outlier values
novec_8(outlierIndex_novec_8) = [];
best_time_novec_8 = max( novec_8 );

percntiles_novec_16 = prctile(novec_16,[5 95]); %5th and 95th percentile
outlierIndex_novec_16 = novec_16 < percntiles_novec_16(1) | novec_16 > percntiles_novec_16(2);
%remove outlier values
novec_16(outlierIndex_novec_16) = [];
best_time_novec_16 = max( novec_16 );

percntiles_novec_32 = prctile(novec_32,[5 95]); %5th and 95th percentile
outlierIndex_novec_32 = novec_32 < percntiles_novec_32(1) | novec_32 > percntiles_novec_32(2);
%remove outlier values
novec_32(outlierIndex_novec_32) = [];
best_time_novec_32 = max( novec_32 );


dataset = [1 2 4 8 16 32];
dataset_novec = [1 2 4 8 16 32];

dataset_pgres_seq = [1 2 4 8 16 32];
dataset_pgres_par = [1 2 4 8 16 32];

time_olap = [ best_time_1 best_time_2 best_time_4 best_time_8 best_time_16 best_time_32 ];
time_olap_novec = [ best_time_novec_1 best_time_novec_2 best_time_novec_4 best_time_novec_8 best_time_novec_16 best_time_novec_32];

time_pgres_seq = [ best_time_pgres_1 best_time_pgres_2 best_time_pgres_4 best_time_pgres_8 best_time_pgres_16 best_time_pgres_32  ];
time_pgres_par = [ best_time_pgres_par_1 best_time_pgres_par_2 best_time_pgres_par_4 best_time_pgres_par_8 best_time_pgres_par_16 best_time_pgres_par_32  ];


loglog(dataset_novec,time_olap_novec,'s--','Color', color9, 'MarkerEdgeColor' , color9, 'MarkerSize', 14);
hold on;

loglog(dataset,time_olap,'s--','Color', color0,'MarkerSize', 14);
hold on;

loglog(dataset_pgres_seq,time_pgres_seq,'d--','Color', color2,'MarkerSize', 14);
hold on;

loglog(dataset_pgres_par,time_pgres_par,'d--','Color', color6,'MarkerSize', 14);
hold on;




B = [ 1 best_time_novec_1 ; 2 best_time_novec_2 ; 4 best_time_novec_4 ; 8 best_time_novec_8 ; 16 best_time_novec_16 ; 32 best_time_novec_32 ];
str = num2str([ best_time_novec_1 ; best_time_novec_2 ; best_time_novec_4 ; best_time_novec_8  ; best_time_novec_16 ; best_time_novec_32 ]);
text(B(:,1)*1.10, B(:,2)*1.10,str);

A = [ 1 best_time_1 ; 2 best_time_2 ; 4 best_time_4 ; 8 best_time_8 ; 16 best_time_16 ; 32 best_time_32 ];
str = num2str([ best_time_1 ; best_time_2 ; best_time_4 ; best_time_8  ; best_time_16 ; best_time_32 ]);
text(A(:,1)*1.10, A(:,2)*1,str);

C = [ 1 best_time_pgres_1 ; 2 best_time_pgres_2 ; 4 best_time_pgres_4 ; 8 best_time_pgres_8  ; 16 best_time_pgres_16  ; 32 best_time_pgres_32  ];
str1 = num2str([ best_time_pgres_1 ; best_time_pgres_2 ; best_time_pgres_4 ; best_time_pgres_8 ; best_time_pgres_16 ; best_time_pgres_32 ]);
text(C(:,1)*1.10, C(:,2)*1.00,str1);

D = [ 1 best_time_pgres_par_1 ; 2 best_time_pgres_par_2 ; 4 best_time_pgres_par_4 ; 8 best_time_pgres_par_8  ; 16 best_time_pgres_par_16  ; 32 best_time_pgres_par_32  ];
str2 = num2str([ best_time_pgres_par_1 ; best_time_pgres_par_2 ; best_time_pgres_par_4 ; best_time_pgres_par_8 ; best_time_pgres_par_16 ; best_time_pgres_par_32 ]);
text(D(:,1)*1.10, D(:,2)*1.10,str2);


grid on;
set(gca, 'YTick', [1 2 4 8 16 32 64 128 ]);

set(gca, 'XTick', [1 2 4 8 16 32 ]);
xlim([1,32]) ;
ylim([0,164]) ;


set(gca,'YTickLabel',num2str(get(gca,'YTick').'));


l = legend( 'Non Vectorized Sequential Linear Algebra Approach', 'Vectorized Sequential Linear Algebra Approach','Sequential PostgreSQL', 'Parallel PostgreSQL (max\_parallel\_degree = 20)' );


  




set(l,'FontSize',12);
ylabel('Time in seconds');

xlabel('TPC-H scale factor');
t = title({'TPC-H benchmark simplified querie-1 time for solution analysis','for different scale factors, between Linear Algebra approach vs Relational Algebra approach,', 'performing a K-Best (K=3,N=50) time measurement technique'},'interpreter','latex')

set(t,'FontSize',24);
set(gca,'fontsize',12);




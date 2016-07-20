
red = [ 255 0 0 ];
green = [ 0 102 0 ];
blue = [ 0 0 255 ];
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
  set(FigHandle, 'Position', [0, 0, 760, 760]);
cd ('../../../src');

csv_results_1 = readtable('timing/timings_vec_nomkl_1.csv','ReadVariableNames',true);
csv_results_2 = readtable('timing/timings_vec_nomkl_2.csv','ReadVariableNames',true);
csv_results_4 = readtable('timing/timings_vec_nomkl_4.csv','ReadVariableNames',true);
csv_results_8 = readtable('timing/timings_vec_nomkl_8.csv','ReadVariableNames',true);
csv_results_16 = readtable('timing/timings_vec_nomkl_16.csv','ReadVariableNames',true);
csv_results_32 = readtable('timing/timings_vec_nomkl_32.csv','ReadVariableNames',true);

time_1 = csv_results_1.total;
time_2 = csv_results_2.total;
time_4 = csv_results_4.total;
time_8 = csv_results_8.total;
time_16 = csv_results_16.total;
time_32 = csv_results_32.total;

best_time_1 = min( time_1 );
best_time_2 = min( time_2 );
best_time_4 = min( time_4 );
best_time_8 = min( time_8 );
best_time_16 = min( time_16 );
best_time_32 = min( time_32 );


best_time_seq_1 = 1.55564;
best_time_seq_2 = 3.22935;
best_time_seq_4= 5.65117;
best_time_seq_8 = 10.5439;
best_time_seq_16 = 21.2006;
best_time_seq_32 = 48.8512;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


csv_results_1 = readtable('timing/timings_test_1.csv','ReadVariableNames',true);
csv_results_2 = readtable('timing/timings_test_2.csv','ReadVariableNames',true);
csv_results_4 = readtable('timing/timings_test_4.csv','ReadVariableNames',true);
csv_results_8 = readtable('timing/timings_test_8.csv','ReadVariableNames',true);
csv_results_16 = readtable('timing/timings_test_16.csv','ReadVariableNames',true);
csv_results_32 = readtable('timing/timings_test_32.csv','ReadVariableNames',true);

time_test_1 = csv_results_1.total;
time_test_2 = csv_results_2.total;
time_test_4 = csv_results_4.total;
time_test_8 = csv_results_8.total;
time_test_16 = csv_results_16.total;
time_test_32 = csv_results_32.total;

best_time_test_1 = min( time_test_1 );
best_time_test_2 = min( time_test_2 );
best_time_test_4 = min( time_test_4 );
best_time_test_8 = min( time_test_8 );
best_time_test_16 = min( time_test_16 );
best_time_test_32 = min( time_test_32 );


%%%%%%%%%%%%%%%%%%%

seq_pgres_1 = csvread('timing/pgres_seq_1.csv') / 1000;
seq_pgres_2 = csvread('timing/pgres_seq_2.csv') / 1000;
seq_pgres_4 = csvread('timing/pgres_seq_4.csv') / 1000;
seq_pgres_8 = csvread('timing/pgres_seq_8.csv') / 1000;
seq_pgres_16 = csvread('timing/pgres_seq_16.csv') / 1000;
seq_pgres_32 = csvread('timing/pgres_seq_32.csv') / 1000;
best_time_pgres_1 = min( seq_pgres_1 );
best_time_pgres_2 = min( seq_pgres_2 );
best_time_pgres_4 = min( seq_pgres_4 );
best_time_pgres_8 = min( seq_pgres_8 );
best_time_pgres_16 = min( seq_pgres_16 );
best_time_pgres_32 = min( seq_pgres_32 );

par_pgres_1 = csvread('timing/pgres_par_1.csv') / 1000;
par_pgres_2 = csvread('timing/pgres_par_2.csv') / 1000;
par_pgres_4 = csvread('timing/pgres_par_4.csv') / 1000;
par_pgres_8 = csvread('timing/pgres_par_8.csv') / 1000;
par_pgres_16 = csvread('timing/pgres_par_16.csv') / 1000;
par_pgres_32 = csvread('timing/pgres_par_32.csv') / 1000;
best_time_pgres_par_1 = min( par_pgres_1 )
best_time_pgres_par_2 = min( par_pgres_2 )
best_time_pgres_par_4 = min( par_pgres_4 )
best_time_pgres_par_8 = min( par_pgres_8 )
best_time_pgres_par_16 = min( par_pgres_16 )
best_time_pgres_par_32 = min( par_pgres_32 )



dataset = [1 2 4 8 16 32];

%time_olap = [ best_time_1 best_time_2 best_time_4 best_time_8 best_time_16 best_time_32 ];

time_pgres_seq = [ best_time_pgres_1 best_time_pgres_2 best_time_pgres_4 best_time_pgres_8 best_time_pgres_16 best_time_pgres_32  ];
time_pgres_par = [ best_time_pgres_par_1 best_time_pgres_par_2 best_time_pgres_par_4 best_time_pgres_par_8 best_time_pgres_par_16 best_time_pgres_par_32  ];

time_olap_seq = [ best_time_seq_1 best_time_seq_2 best_time_seq_4 best_time_seq_8 best_time_seq_16 best_time_seq_32 ];
time_olap_new = [ best_time_1 best_time_2 best_time_4 best_time_8 best_time_16 best_time_32 ];
time_olap_new_v1 = [ best_time_test_1 best_time_test_2 best_time_test_4 best_time_test_8 best_time_test_16 best_time_test_32 ];


%loglog(dataset,time_olap,'s--','Color', color0,'MarkerSize', 14);
%hold on;

loglog(dataset,time_olap_seq,'s--','Color', color0,'LineWidth',2 , 'MarkerSize', 14);
hold on;

loglog(dataset,time_olap_new,'o--','Color', color7, 'LineWidth',2 , 'MarkerSize', 14);
hold on;

loglog(dataset,time_olap_new_v1,'+--','Color', color9, 'LineWidth',2 , 'MarkerSize', 14);
hold on;

loglog(dataset,time_pgres_seq,'d--','Color', color2, 'LineWidth',2 ,'MarkerSize', 14);
hold on;

loglog(dataset,time_pgres_par,'x--','Color', color6,'LineWidth',2 ,'MarkerSize', 14);
hold on;

A = [ 1 best_time_1 ; 2 best_time_2 ; 4 best_time_4 ; 8 best_time_8 ; 16 best_time_16 ; 32 best_time_32 ];
str = num2str([ best_time_1 ; best_time_2 ; best_time_4 ; best_time_8  ; best_time_16 ; best_time_32 ]);
text(A(:,1)*1.10, A(:,2)*1,str);

B = [ 1 best_time_test_1 ; 2 best_time_test_2 ; 4 best_time_test_4 ; 8 best_time_test_8 ; 16 best_time_test_16 ; 32 best_time_test_32 ];
str = num2str([ best_time_test_1 ; best_time_test_2 ; best_time_test_4 ; best_time_test_8  ; best_time_test_16 ; best_time_test_32 ]);
text(B(:,1)*1.10, B(:,2)*1,str);

C = [ 1 best_time_pgres_1 ; 2 best_time_pgres_2 ; 4 best_time_pgres_4 ; 8 best_time_pgres_8  ; 16 best_time_pgres_16  ; 32 best_time_pgres_32  ];
str1 = num2str([ best_time_pgres_1 ; best_time_pgres_2 ; best_time_pgres_4 ; best_time_pgres_8 ; best_time_pgres_16 ; best_time_pgres_32 ]);
text(C(:,1)*1.10, C(:,2)*1.00,str1);

D = [ 1 best_time_pgres_par_1 ; 2 best_time_pgres_par_2 ; 4 best_time_pgres_par_4 ; 8 best_time_pgres_par_8  ; 16 best_time_pgres_par_16  ; 32 best_time_pgres_par_32  ];
str2 = num2str([ best_time_pgres_par_1 ; best_time_pgres_par_2 ; best_time_pgres_par_4 ; best_time_pgres_par_8 ; best_time_pgres_par_16 ; best_time_pgres_par_32 ]);
text(D(:,1)*1.10, D(:,2)*1.10,str2);

Z = [ 1 best_time_seq_1 ; 2 best_time_seq_2 ; 4 best_time_seq_4 ; 8 best_time_seq_8 ; 16 best_time_seq_16 ; 32 best_time_seq_32 ];
str = num2str([ best_time_seq_1 ; best_time_seq_2 ; best_time_seq_4 ; best_time_seq_8  ; best_time_seq_16 ; best_time_seq_32 ]);
text(Z(:,1)*1.10, Z(:,2)*1,str);

grid on;
set(gca, 'YTick', [1 2 4 8 16 32 64 128 ]);

set(gca, 'XTick', [1 2 4 8 16 32 ]);
xlim([1,32]) ;
ylim([0,64]) ;


set(gca,'YTickLabel',num2str(get(gca,'YTick').'));


l = legend( 'Vectorized Sequential Linear Algebra Approach',  'Parallel Linear Algebra Approach' , 'Parallel Linear Algebra Approach New Version',  'Sequential PostgreSQL', 'Parallel PostgreSQL (max\_parallel\_degree = 20)' );


set(l,'FontSize',12);
ylabel('Time in seconds');

xlabel('TPC-H scale factor');
t = title({'TPC-H benchmark simplified query-1 time for solution analysis','for different scale factors, between Linear Algebra approach vs Relational Algebra approach,', 'performing a K-Best (K=3,N=50) time measurement technique'},'interpreter','latex');
set(gca,'fontsize',12);

set(t,'FontSize',14);





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


speedup_olap = time_olap_seq ./ time_olap_new_v1;

speedup_vs_pgres = time_pgres_par ./ time_olap_new_v1;

%loglog(dataset,time_olap,'s--','Color', color0,'MarkerSize', 14);
%hold on;


loglog(dataset,speedup_olap,'o--','Color', color7, 'LineWidth',2 , 'MarkerSize', 14);
hold on;

loglog(dataset,speedup_vs_pgres,'x--','Color', color6,'LineWidth',2 ,'MarkerSize', 14);
hold on;


grid on;
set(gca, 'YTick', [1 2 4 8 16 32 64 128 ]);
set(gca, 'XTick', [1 2 4 8 16 32 ]);
xlim([1,32]) ;
ylim([0,24]) ;

set(gca,'YTickLabel',num2str(get(gca,'YTick').'));

l = legend( 'Speedup SLA vs PLA',  'Speedup PLA vs PRA'  );

set(l,'FontSize',12);
ylabel('Speedup');

xlabel('TPC-H scale factor');
t = title({'TPC-H benchmark simplified query-1 speedup analysis','for different scale factors, between Sequential Linear Algebra vs Parallel Linear Algebra,','and Parallel Linear Algebra vs Parallel Relational Algebra'},'interpreter','latex');
set(gca,'fontsize',12);

set(t,'FontSize',14);




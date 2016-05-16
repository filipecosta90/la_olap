  FigHandle = figure;
  set(FigHandle, 'Position', [0, 0, 640, 480]);
cd ('../../../src/timing');
dataset_1 = csvread('timings_search_1.dat');
dataset_2 = csvread('timings_search_2.dat');
dataset_4 = csvread('timings_search_4.dat');
dataset_8 = csvread('timings_search_8.dat');
dataset_16 = csvread('timings_search_16.dat');


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

dataset = [1 2 4 8 16];
time_olap = [ best_time_1 best_time_2 best_time_4 best_time_8 best_time_16 ];


loglog(dataset,time_olap,'ro--','MarkerSize', 12);
hold on;
%loglog(threads_cg_c_641_93,tempo_cg_c_641_93,'r+--','Color', cores(2,:),'MarkerSize', 12);
%hold on;

%%%%%


grid on;
set(gca, 'XTick', [1 2 4 8 16 32 ]);
xlim([1,32]) ;
%ylim([0,1000]) ;


set(gca,'YTickLabel',num2str(get(gca,'YTick').'));


l = legend('431 - gcc 4.9.0 -O3','641 - gcc 4.9.0 -O3' , '652 - gcc 4.9.0 -O3', '662 - gcc 4.9.0 -O3');


  




set(l,'FontSize',12);
ylabel('Tempo (segundos)');

xlabel('Num. Threads OpenMP');
t = title({'Rela\c{c}\~ao entre Tempo Total em segundos para o kernel OMP - CG','Classe de dados C para compilador gcc 4.9.0 com flags de compila\c{c}\~ao -O3'},'interpreter','latex')

set(t,'FontSize',24);
set(gca,'fontsize',12);




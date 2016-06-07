 

cpu_csv = readtable('__cpu_usage/CPU_postgressql_dstat.csv','ReadVariableNames',false);
disk_csv = readtable('__disk_usage/DISK_postgressql_dstat.csv','ReadVariableNames',false);
memory_csv = readtable('__memory_usage/MEMORY_postgressql_dstat.csv','ReadVariableNames',false);
net_csv = readtable('__net_usage/NET_postgressql_dstat.csv','ReadVariableNames',false);
system_csv = readtable('__system_usage/SYSTEM_postgressql_dstat.csv','ReadVariableNames',false);

%%%%%% CPU STATS %%%%%
cpu_usr = table2array( cpu_csv ( :, [3])); 
cpu_usr = cpu_usr (1:end-4,:)

cpu_sys = table2array( cpu_csv ( :, [4])); 
cpu_sys = cpu_sys (1:end-4,:)

cpu_idl = table2array( cpu_csv ( :, [5])); 
cpu_idl = cpu_idl (1:end-4,:)

cpu_wait = table2array( cpu_csv ( :, [6]));
cpu_wait = cpu_wait (1:end-4,:)


%%%%%% DISK STATS %%%%%
disk_read = table2array( disk_csv ( :, [3])); 
disk_read = disk_read (1:end-4,:)

disk_write = table2array( disk_csv ( :, [4]));
disk_write = disk_write (1:end-4,:)


%%%%%% MEMORY STATS %%%%%
memory_used = table2array( memory_csv ( :, [3])); 
memory_used = memory_used (1:end-4,:)

memory_free = table2array( memory_csv ( :, [4])); 
memory_free = memory_free (1:end-4,:)

memory_total = memory_used + memory_free;
memory_used = memory_used / memory_total * 100;
memory_free = memory_free / memory_total * 100;


%%%%%% NET STATS %%%%%
net_recv = table2array( net_csv ( :, [3]) ); 
net_recv = net_recv (1:end-4,:)

net_send = table2array( net_csv ( :, [4]) ); 
net_send = net_send (1:end-4,:)

%%%%%% SYSTEM STATS %%%%%
system_int = table2array( system_csv ( :, [3])); 
system_int = system_int (1:end-4,:)

system_csw = table2array( system_csv ( :, [4])) ; 
system_csw = system_csw (1:end-4,:)


bg = [1 1 1; 0 0 0];
cores = distinguishable_colors(100,bg);

 figure (1)


 hFig = figure(1);
 set(gcf,'PaperPositionMode','auto')
 set(hFig, 'Position', [0 0 960 960])
set(gca,'Unit','normalized','Position',[0 0 1 1]);

subplot(3,2,1)       % add first plot in 2 x 2 grid
plot(system_int);
hold on;
plot(system_csw);
hold on;
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));
l = legend('System INT','System CSW' );
ylabel('N. de ocorr\^encias', 'interpreter','latex');
xlabel('Tempo (s)');

title({'System Stats'},'interpreter','latex')


subplot(3,2,2)       % add second plot in 2 x 2 grid
plot(disk_read);
hold on;
plot(disk_write);
hold on;
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));
l = legend('Disk Read','Disk Write' );
ylabel('Bytes');
xlabel('Tempo (s)');

title({'Disk Stats'},'interpreter','latex')

subplot(3,2,3)       % add third plot in 2 x 2 grid

plot(net_recv);
hold on;
plot(net_send);
hold on;
set(gca,'YTickLabel',num2str(get(gca,'YTick').'));
l = legend('Net Receive','Net Send' );
ylabel('Bytes');
xlabel('Tempo (s)');

title({'Net Stats'},'interpreter','latex')


subplot(3,2,4)       % add fourth plot in 2 x 2 grid

mem = [ memory_used memory_free ];
h = bar(mem,'stacked');
set(h(1),'facecolor','r');
set(h(2),'facecolor','g');



set(gca,'YTickLabel',num2str(get(gca,'YTick').'));
l = legend('% Memory Used','% Memory Free' );
ylabel('Percentagem');
xlabel('Tempo (s)');
ylim([0 100]) ;


title({'Memory Stats'},'interpreter','latex');

subplot(3,2,[5 6])    % add third plot to span positions 3 and 4

cpu = [ cpu_usr cpu_sys cpu_idl cpu_wait ];

hb = bar(cpu,'stacked');
ylim([0 100]) ;
title({'CPU Stats'},'interpreter','latex');
xlabel('Tempo (s)');
ylabel('Percentagem');
cores = distinguishable_colors(8,bg);
set(hb(1), 'FaceColor',cores(1,:));
set(hb(2), 'FaceColor',cores(3,:));
set(hb(3), 'FaceColor',cores(8,:));
set(hb(4), 'FaceColor',cores(2,:));
l = legend('% CPU USR','% CPU SYS','% CPU IDL','% CPU WAIT' );








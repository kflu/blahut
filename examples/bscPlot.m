ce0 = load('ce0.txt');
ce1 = load('ce0.1.txt');
ce2 = load('ce0.2.txt');
ce3 = load('ce0.3.txt');
ce4 = load('ce0.4.txt');
ce5 = load('ce0.5.txt');

figure(99);
clg;
hold on;
grid on;

title('Capacity Expense Curve of Binary Symmetric Channel');
xlabel('Expense');
ylabel('Capacity (bits/c.c)')

plot(ce0(:,1),ce0(:,2), ';BSC, p=0.0;r-');
plot(ce1(:,1),ce1(:,2), ';BSC, p=0.1;g-');
plot(ce2(:,1),ce2(:,2), ';BSC, p=0.2;b-');
plot(ce3(:,1),ce3(:,2), ';BSC, p=0.3;k-');
plot(ce4(:,1),ce4(:,2), ';BSC, p=0.4;c-');
plot(ce5(:,1),ce5(:,2), ';BSC, p=0.5;m-');

pause;

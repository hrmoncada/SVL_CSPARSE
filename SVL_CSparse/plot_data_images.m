%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all   % close all open such as : figures, fuctions, etc
clc         % clear the command prompt
clear all   % clear all variables

% loading data 
load data1
load data2
load data3
load data4
load data5
load data6
load data7
load data8
load data9
load data10
%  load data11
%  load data12
%  load data13
%  load data14
load data11
load data12
load data13

M1 = length(data1(:,1)); %  U1 zeros matrix      (Nx x Ny)
M2 = length(data2(:,1)); %  U2 + Triangle        (Nx x Ny)
M3 = length(data3(:,1)); %  A_re = real(FFT(U))  (Nx x Ny)
M4 = length(data4(:,1)); %  A_im = img(FFT(U))   (Nx x Ny)
M5 = length(data5(:,1)); %  SWAP(A_re)           (Nx x Ny)
M6 = length(data6(:,1)); %  SWAP(A_im)           (Nx x Ny)
M7 = length(data7(:,1)); %  real(TA)             (NM x NN)
M8 = length(data8(:,1)); %  imag(TA)             (NM x NN)
M9 = length(data9(:,1)); %  RQS                  (10*NPx x  10*NPy) =   (New_Nx x New_Ny)
M10 = length(data10(:,1)); %  THETA              (New_Nx x New_Ny)
%M11 = length(data11(:,1)); %  DX                (New_Nx x New_Ny)
%M12 = length(data12(:,1)); %  D2X               (New_Nx x New_Ny)
%M13 = length(data13(:,1)); %  DY                (New_Nx x New_Ny)
%M14 = length(data14(:,1)); %  D2Y               (New_Nx x New_Ny)
M15 = length(data11(:,1)); %  PHI                (M x M)
M16 = length(data12(:,1)); %  S                  (M x M)
M17 = length(data13(:,1)); %  UC                 (M x M)

%  
for j = 1 : M1     U1(:,j) =  data1(:,j); end
for j = 1 : M2     U2(:,j) =  data2(:,j); end	
for j = 1 : M3   A_re(:,j) =  data3(:,j); end
for j = 1 : M4   A_im(:,j) =  data4(:,j); end
for j = 1 : M5  SA_re(:,j) =  data5(:,j); end
for j = 1 : M6  SA_im(:,j) =  data6(:,j); end
for j = 1 : M7  TA_re(:,j) =  data7(:,j); end
for j = 1 : M8  TA_im(:,j) =  data8(:,j); end
for j = 1 : M9    RSQ(:,j) =  data9(:,j); end
for j = 1 : M10 THETA(:,j) =  data10(:,j); end
%  for j = 1 : M11    DX(:,j) =  data10(:,j); end
%  for j = 1 : M12   D2X(:,j) =  data11(:,j); end
%  for j = 1 : M13    DY(:,j) =  data12(:,j); end
%  for j = 1 : M14   D2Y(:,j) =  data13(:,j); end
for j = 1 : M15    PHI(:,j) =  data11(:,j); end
for j = 1 : M16      S(:,j) =  data12(:,j); end
for j = 1 : M17     UC(:,j) =  data13(:,j); end

% Grid size
Lx = 1; % x-axis unit cell grid size
Ly = 1; % y-axis unit cell grid size


% Step size

dx = Lx/(M1-1); % x-axis step size
dy = Ly/(M1-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(1)
imagesc(x, y, U1)      %U1 zeros matrix (Nx x Ny)
axis equal tight
%grid on
%print (1,"Fig_1_Empty_unit_cell.epsc");
print -depsc Fig_1_Empty_unit_cell.eps;

%colormap(hot)
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('Zeros array')

% Step size
dx = Lx/(M2-1); % x-axis step size
dy = Ly/(M2-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
U2 = rot90(U2,-1);
figure(2)
imagesc(x, y, U2)          %U2 + Triangle(Nx x Ny)
axis equal tight
%grid on
%print (2,"Fig_1_empty_unit_cell.epsc");
print -depsc Fig_2_Device_unit_cell.eps;

%colormap(hot)
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('Triangle')

% Step size
dx = Lx/(M3-1); % x-axis step size
dy = Ly/(M3-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(3)
imagesc(x, y, A_re)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('FFT real ')
%print (3,"Fig_3_FFT_Real.epsc");
print -depsc Fig_3_FFT_Real.eps;

% Step size
dx = Lx/(M4-1); % x-axis step size
dy = Ly/(M4-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(4)
imagesc(x, y, A_im)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('FFT imaginary ')
%print (4,"Fig_4_FFT_Imag.epsc");
print -depsc Fig_4_FFT_Imag.eps;

% Step size
dx = Lx/(M5-1); % x-axis step size
dy = Ly/(M5-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(5)
imagesc(x, y, SA_re)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('SWAP(FFT(A_re) ')
%print (5,"Fig_5_Swap_Real.epsc");
print -depsc Fig_5_Swap_Real.eps;


% Step size
dx = Lx/(M6-1); % x-axis step size
dy = Ly/(M6-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(6)
imagesc(x, y, SA_im)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('SWAP(FFT(A_im) ')
%print (6,"Fig_6_Swap_Imag.epsc");
print -depsc Fig_6_Swap_Imag.eps;

% Step size
dx = Lx/(M7-1); % x-axis step size
dy = Ly/(M7-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(7)
imagesc(x, y, TA_re)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('AT_re')
%print (7,"Fig_7_Trucate_Real.epsc");
print -depsc Fig_7_Trucate_Real.eps;


% Step size
dx = Lx/(M8-1); % x-axis step size
dy = Ly/(M8-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(8)
imagesc(x, y, TA_im)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('AT_im')
%print (8,"Fig_8_Trucate_Imag.epsc");
print -depsc Fig_8_Trucate_Imag.eps;

dx = Lx/(M9-1); % x-axis step size
dy = Ly/(M9-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(9)
imagesc(RSQ)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('RSQ')
%print (9,"Fig_9_RSQ.epsc");
print -depsc Fig_9_RSQ.eps;

dx = Lx/(M10-1); % x-axis step size
dy = Ly/(M10-1); % x-axis step size
x = 0 : dx : 1;
y = 0 : dx : 1;
figure(10)
imagesc(THETA)
axis equal tight
%grid on
% axis manual
% axis fill
%colorbar
%xlabel('x')
%ylabel('y')
%title('THETA')
%print (10,"Fig_10_Theta.epsc");
print -depsc Fig_10_Theta.eps;

%  figure(2)
%  subplot(2,2,1);
%  imagesc(DX)
%  axis equal tight
%  grid on
%  title('DX')
%  
%  subplot(2,2,2);
%  imagesc(DX)
%  axis equal tight
%  grid on
%  title('D2X')
%  
%  subplot(2,2,3);
%  imagesc(DY)
%  axis equal tight
%  grid on
%  title('DY')
%  
%  subplot(2,2,4);
%  imagesc(D2Y)
%  axis equal tight
%  grid on
%  title('D2Y')

  figure(11)
  dx = Lx/(M15-1); % x-axis step size
  dy = Ly/(M15-1); % x-axis step size
  x = 0 : dx : 1;
  y = 0 : dx : 1;
  subplot(1,3,1);
  imagesc(PHI)
  axis equal tight
  grid on
  title('PHI')

  dx = Lx/(M16-1); % x-axis step size
  dy = Ly/(M16-1); % x-axis step size
  x = 0 : dx : 1;
  y = 0 : dx : 1; 
  subplot(1,3,2);
  imagesc(S)
  axis equal tight
  grid on
  title('S')

  dx = Lx/(M17-1); % x-axis step size
  dy = Ly/(M17-1); % x-axis step size
  x = 0 : dx : 1;
  y = 0 : dx : 1;
  subplot(1,3,3);
  imagesc(UC)
  axis equal tight
  grid on
  title('UC')
  

%print (11,"Figure_11.eps");
print -depsc Figure_11.eps;




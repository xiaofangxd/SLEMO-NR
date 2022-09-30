cd(fileparts(mfilename('fullpath')));
addpath(genpath(cd));
clear;
clc;
pop_all_4 = load('pop_all_4.mat');
pop_all_5 = load('pop_all_5.mat');
pop_all_7 = load('pop_all_7.mat');
fit_all_6 = load('fit_all_6.mat');
fit_all_1 = load('fit_all_1.mat');
fit_all_2 = load('fit_all_2.mat');
fit_all_3 = load('fit_all_3.mat');


choose = 9;

% CS 1
f11 = fit_all_1.fit_all(choose,1);
f12 = fit_all_1.fit_all(choose,2);

% LASSO 2
f21 = fit_all_2.fit_all(choose,1);
f22 = fit_all_2.fit_all(choose,2);

% OMP 3
f31 = fit_all_3.fit_all(choose,1);
f32 = fit_all_3.fit_all(choose,2);

% SparseEA 4
population = pop_all_4.pop_all(choose,:);
NonDominated = NDSort(population.objs,1) == 1;
NondomPopulation = population(NonDominated);
f4 = NondomPopulation.objs;
f41 = f4(:,1);
f42 = f4(:,2);

% MOEAPSL 5
population = pop_all_5.pop_all(choose,:);
NonDominated = NDSort(population.objs,1) == 1;
NondomPopulation = population(NonDominated);
f5 = NondomPopulation.objs;
f51 = f5(:,1);
f52 = f5(:,2);

% MAST-Net 6
f61 = fit_all_6.fit_all(choose,1);
f62 = fit_all_6.fit_all(choose,2);

% SLEMO-NR 7
population = pop_all_7.pop_all(choose,:);
NonDominated = NDSort(population.objs,1) == 1;
NondomPopulation = population(NonDominated);
f7 = NondomPopulation.objs;
f71 = f7(:,1);
f72 = f7(:,2);



% f1 = figure;
plot(f12,f11,'d',f22,f21,'s',f32,f31,'+',f42,f41,'^',f52,f51,'o',f62,f61,'p',f72,f71,'*');

xlabel('||\bf{\it{X_i}}\rm||_0','interpreter','tex');
ylabel('\bf{\it{M}}\rm^{-1}||\bf{\it{A_iX_i}}\rm{-}\bf{\it{Y_i}}\rm||^2_2','Interpreter','tex');
legend('CS','LASSO','OMP','SparseEA','MOEA/PSL','MAST-Net','SLEMO-NR');
set(gca,'Fontname','times new Roman','Fontsize',12);
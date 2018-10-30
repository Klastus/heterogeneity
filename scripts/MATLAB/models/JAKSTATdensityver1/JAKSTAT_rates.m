function R = JAKSTAT_rates(x, par, t, stimulus)
R =[
%(par(1)*x(15))*par(5)*stimulus(1)/(1+par(5)*stimulus(1))*x(1)/(par(6)+x(1));
par(1)*(par(13)*x(16)+par(13)*x(17))*par(11)*x(1)/(par(11)*par(6)+par(11)*x(1));
par(2)*par(11)*x(2)*par(11)*x(2);
par(3)*par(11)*x(3);
par(4)*par(12)*x(4);
par(4)*par(12)*x(5);
par(4)*par(12)*x(6);
par(4)*par(12)*x(7);
par(4)*par(12)*x(8);
par(4)*par(12)*x(9);
par(4)*par(12)*x(10);
par(4)*par(12)*x(11);
par(4)*par(12)*x(12);
par(4)*par(12)*x(13);
par(5)*par(7)*stimulus(1)*par(13)*x(15);
par(8)*par(13)*x(16);
par(9)*par(13)*x(16);
par(10)*par(13)*x(17);
];
 end

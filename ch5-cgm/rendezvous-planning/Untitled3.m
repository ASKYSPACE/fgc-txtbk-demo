
subplot(1,4,1)
h1=openfig('T200obsoff_T10.fig','reuse');
h2=openfig('T200obsoff_T42.fig','reuse');
h3=openfig('T200obsoff_T82.fig','reuse');
h4=openfig('T200obsoff_T82.fig','reuse');
figure(12)
s1=subplot(1,4,1);
copyobj(get(get(h1,'Children'),'Children'),s1)
s2=subplot(1,4,2);
copyobj(get(get(h2,'Children'),'Children'),s2)
s3=subplot(1,4,3);
copyobj(get(get(h3,'Children'),'Children'),s3)
s4=subplot(1,4,4);
copyobj(get(get(h4,'Children'),'Children'),s4)
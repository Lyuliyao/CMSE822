N = 4;
x = (0:N)/N
for i = 1:N-1
    y = (x==i/N);
    plot(x,y,'DisplayName',['$\varphi$',num2str(i)])
    hold on
end
hold off
legend('Interpreter','latex');
print basis_figure.eps -depsc2 -r600;
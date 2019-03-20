function KapitzaPendulum(;l=1, g=9.8, ax=0, ay=0, omegax=0, omegay=0, phix=0, phiy=0)
    """
    A Kapitza Pendulum function defined by the parameters.
    
    :param l: length.
    :param g: acceleration due to gravity.
    
    :params ax, omegax, phix: the x coordinate of the pivot xp(t) = ax cos(omegax t + phix).
    :params ay, omegay, phiy: the y coordinate of the pivot yp(t) = ay cos(omegay t + phiy).
    
    :return: the function f where \vec{theta}'(t) = f(t, \vec{theta(t)}), x, y.
    """
    function f(t, theta)
        """
        For [theta'(t), theta''(t)] = f(t, [theta(t), theta'(t)]).
        
        :param t: time.
        :param theta: vector [theta(t), theta'(t)].
        
        :return: \vec{theta'(t)}.
        """
        th, thdot = theta
        [thdot, ax/l*omegax^2*cos(omegax*t+phix)*cos(th)-(g-ay*omegay^2*cos(omegay*t+phiy))*sin(th)/l]
    end
    function x(t, theta)
        """ x(t) = xp(t) + l sin(theta(t)) """
        ax*cos(omegax*t+phix)+l*sin(theta)
    end
    function y(t, theta)
        """ y(t) = yp(t) - y cos(theta(t)) """
        ay*cos(omegay*t+phiy)-l*cos(theta)
    end
    
    f, x, y
end

#Default arguments represents simple pendulum.
gfunc, x, y = KapitzaPendulum(ax=.2, omegax=100, ay=.2, omegay=100, phiy=pi/2.0)
t, tmax, h = 0.0, 2, 0.0001
y0 = [pi/4.0, 0.0]
ts, theta = MultiStepper(gfunc, AdamsMoult4, t, y0, tmax, h)

using PyPlot

# Theta plot
figure(0)
plot(ts, [theta[i][1] for i=1:length(theta)], "o-", label="Position")
plot(ts, [theta[i][2] for i=1:length(theta)], "x--", label="Velocity")
xlabel("time t")
ylabel("angle theta (radians)")
title("Kapitza Pendulum")
legend()

# Parametric plot
figure(1)
xaxis = [x(ts[i], theta[i][1]) for i=1:length(ts)]
yaxis = [y(ts[i], theta[i][1]) for i=1:length(ts)]
plot(xaxis, yaxis, "o-", label="Parametric Plot")
title("Kapitza Pendulum")
xlabel("x")
ylabel("y")
legend()
#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QDesktopWidget>
#include <math.h>

double f_V_x( double V_x, double k, double m)
{
	return -k/m*V_x;
}

double f_V_y( double V_y, double k, double m, double g)
{
	return -k/m*V_y- g;
}

double f_Vx1( double V_x1, double V_y1, double k, double m)
{
	return -k*V_x1*sqrt(V_x1*V_x1 + V_y1*V_y1)/m;
}

double f_Vy1( double V_x1, double V_y1, double k, double m, double g)
{
	return -k*V_y1*sqrt(V_x1*V_x1 + V_y1*V_y1)/m - g;
}


std::pair<double, double> position_in_space (double m, double k, double g, double t, double V_x0, double V_y0)
{
	double V_x = V_x0*exp(-k/m * t),
		   V_y = V_y0*exp(-k/m * t) - g*m/k*(1 - exp(-k/m* t));
	std::pair<double, double> coord;
	coord.first = m*V_x/k*(1-exp(-k*t/m)); //x0
	coord.second = m/k*((V_y+m*g/k)*(1-exp(-k*t/m)) - g*t);	//y0

	return coord;
}


QScatterSeries* MainWindow::create_ser_without_resistance( double V_x, double V_y, double g)
{
	QScatterSeries* series = new QScatterSeries();
	series->setMarkerSize(10);
	series->setName("без учета сопротивления воздуха");

	for (double t = 0; t < 20; t += 0.01)
		series->append( V_x*t,V_y*t-g*t*t/2);

	return series;
}


//QLineSeries* MainWindow::create_exact_solution ( double m, double V_x, double V_y, double k, double g)
//{
//	QLineSeries* series = new QLineSeries();
//	series->setName("с учетом сопротивления(точное решение)");

//	for (double t = 0; t < 20; t += 0.01)
//		series->append( m*V_x/k*(1-exp(-k*t/m)), m/k*((V_y+m*g/k)*(1-exp(-k*t/m)) - g*t));

//	return series;
//}


QLineSeries* MainWindow::create_numer_approx ( const int N, double _h, double k, double m, double g, double V_x, double V_y)
{
	QLineSeries* series = new QLineSeries();
	series->setName("с учетом сопротивления(численное приближение)");

	double k1 = 0, k2 = 0, k3 = 0, k4 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
	int i = 1;
	double *pY = new double [N]; pY[0] = 0;
	double *pX = new double [N]; pX[0] = 0;
	double *pV_x = new double [N]; pV_x[0] = V_x;		//массив - скорость: проекция на горизонталь
	double *pV_y = new double [N]; pV_y[0] = V_y;		//массив - скорость: проекция на вертикаль
	double *pV = new double [N]; pV[0] = sqrt(pV_x[0]*pV_x[0] + pV_y[0]*pV_y[0]);

	while (pY[i-1]>-0.00001 && i < N)
		{
			k1 = _h * f_V_x( pV_x[i-1],k, m);
			l1 = _h * f_V_y( pV_y[i-1], k, m, g);
			k2 = _h * f_V_x( pV_x[i-1]+k1/2, k, m);
			l2 = _h * f_V_y( pV_y[i-1]+l1/2, k, m, g);
			k3 = _h * f_V_x( pV_x[i-1]+k2/2, k, m);
			l3 = _h * f_V_y( pV_y[i-1]+l2/2, k, m, g);
			k4 = _h * f_V_x( pV_x[i-1]+k3, k, m);
			l4 = _h * f_V_y( pV_y[i-1]+l3, k, m, g);


			// получаем проекции скорости из предварительно полученных приближенных значений
			pV_x[i] = pV_x[i-1] + ( k1+2*k2+2*k3+k4)/6;
			pV_y[i] = pV_y[i-1] + ( l1+2*l2+2*l3+l4)/6;

			pV[i] = sqrt( pV_x[i]*pV_x[i] + pV_y[i]*pV_y[i] ); //величина вектора скорости из проекций

			pY[i] = pY[i-1]+pV_y[i]*_h;
			pX[i] = pX[i-1]+pV_x[i]*_h; //пересчет координат

			series->append(pX[i], pY[i]);

			i++;
		}
	delete [] pY;
	delete [] pX;
	delete [] pV_x;
	delete [] pV_y;
	delete [] pV;

	return series;
}


QLineSeries* MainWindow::create_relativism ( double _h, double V_y, double V_x,
											 double m, double k, double g, int N)
{
	QLineSeries* series = new QLineSeries();
	series->setName("на релятивистских скоростях(численное приближение)");

	double k1 = 0, k2 = 0, k3 = 0, k4 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
	int j = 1;
	double *pY = new double[N]; pY[0] = 0;
	double *pX = new double[N]; pX[0] = 0;
	double *pV_x = new double [N]; pV_x[0] = V_x;	//массив - скорость: проекция на горизонталь
	double *pV_y = new double [N]; pV_y[0] = V_y;	//массив - скорость: проекция на вертикаль
	double *pV = new double[N]; pV[0] = sqrt(pV_x[0]*pV_x[0] + pV_y[0]*pV_y[0]);

	while (pY[j-1]>-0.00001 && j < N )
		{
			k1 = _h * f_Vx1( pV_y[j-1], pV_x[j-1], k, m);
			l1 = _h * f_Vy1( pV_y[j-1], pV_x[j-1], k, m, g);
			k2 = _h * f_Vx1( pV_y[j-1]+l1/2, pV_x[j-1]+k1/2, k, m);
			l2 = _h * f_Vy1( pV_y[j-1]+l1/2, pV_x[j-1]+k1/2, k, m, g);
			k3 = _h * f_Vx1( pV_y[j-1]+l2/2, pV_x[j-1]+k2/2, k, m);
			l3 = _h * f_Vy1( pV_y[j-1]+l2/2, pV_x[j-1]+k2/2, k, m, g);
			k4 = _h * f_Vx1( pV_y[j-1]+l3, pV_x[j-1]+k3, k, m);
			l4 = _h * f_Vy1( pV_y[j-1]+l3, pV_x[j-1]+k3, k, m, g);

			// получаем проекции скорости из предварительно полученных приближенных значений
			pV_x[j] = pV_x[j-1] + (k1+2*k2+2*k3+k4)/6;
			pV_y[j] = pV_y[j-1] + (l1+2*l2+2*l3+l4)/6;

			pV[j] = sqrt( pV_x[j]*pV_x[j] + pV_y[j]*pV_y[j]); //величина вектора скорости из проекций

			pY[j] = pY[j-1]+pV_y[j]*_h;
			pX[j] = pX[j-1]+pV_x[j]*_h; //пересчет координат

			series->append(pX[j], pY[j]);

			j++;
		}

	delete[] pY;
	delete[] pX;
	delete[] pV;
	delete[] pV_x;
	delete[] pV_y;

	return series;
}

//QLineSeries* MainWindow::create_slow_change (double m_p, double k, double V_xp,
//											 double V_yp, double m_t, double g, double V_xt, double V_yt)
//{
//	QLineSeries* series = new QLineSeries;
//	series->setName("случай постепенного изменения массы: ");

//	double dm;



//	for (double t = 0; t < 20; t += 0.01) {
//	series = create_arbitrary_angle( m_p, k, V_xp, V_yp, dm, g, V_xt, V_yt, );
//	//++=++
//	}
//	return series;
//}


//QLineSeries* MainWindow::create_detaching_mass ( double m_p, double k, double V_xp,
//												 double V_yp, double m_p2, double g, double V_t)
//{
//	QLineSeries* series = new QLineSeries();
//	series->setName("случай резкого изменения массы");
//	bool detached = false;
//	double x0 = 0, y0 = 0;
//	double dt = 0;

//	for (double t = 0; t < 20; t += 0.01) {
//		series->append(x0 + m_p*V_xp/k*(1-exp(-k*(t-dt)/m_p)), y0 + m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*(t-dt)/m_p)) - g*(t-dt)));

//		if(t >= 0.6 && !detached){
//			double V_x0 = V_xp*exp(-k/m_p * t);
//			double V_y0 = V_yp*exp(-k/m_p * t) - g*m_p/k*(1 - exp(-k/m_p* t));
//			x0 = m_p*V_xp/k*(1-exp(-k*t/m_p));
//			y0 = m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*t/m_p)) - g*t);

//			V_xp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_x0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
//			V_yp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_y0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
//			m_p = m_p2;
//			detached = true;
//			dt = t;
//		}
//	}

//	return series;
//}


QLineSeries* MainWindow::create_arbitrary_angle (double m_p, double k, double V_xp,
												 double V_yp, double m_p2, double g, double V_xt, double V_yt,
												 double t_detach, double x0, double y0)
{
	QLineSeries* series = new QLineSeries;

	bool detached = false;
	double dt = 0;

	for (double t = 0; t < 20; t += 0.01) {
		series->append(x0 + m_p*V_xp/k*(1 - exp(-k*(t-dt)/m_p)), y0 + m_p/k*((V_yp+m_p*g/k)*(1 - exp(-k*(t-dt)/m_p)) - g*(t-dt)));

		if(t >= t_detach && !detached){
			double V_x0 = V_xp*exp(-k/m_p * t);
			double V_y0 = V_yp*exp(-k/m_p * t) - g*m_p/k*(1 - exp(-k/m_p* t));
			x0 = m_p*V_xp/k*(1-exp(-k*t/m_p));
			y0 = m_p/k*((V_yp+m_p*g/k)*(1 - exp(-k*t/m_p)) - g*t);

			V_xp = (m_p*V_x0 - (m_p - m_p2)*V_xt)/m_p2;
			V_yp = (m_p*V_y0 - (m_p - m_p2)*V_yt)/m_p2;
			m_p = m_p2;
			detached = true;
			dt = t;
		}
	}

	return series;
}

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);
	QChartView *chartView = new QChartView(this);
	QSize t = this->minimumSize();
	chartView->setMinimumSize(t.width() / 2, t.height());

	chrt = new MyChart(0, 13);

	chartView->setChart(chrt);
	chrt->setTitle("Графики движения тела под углом к горизонту:");


	double m, V_x, V_y, C, S, P, k, g;

	P = 1.29;    // (кг/м^3)
	S = 0.0314;  // S = M_PI*r^2
	C = 0.47;    // безразмерный коэффициент сопротивлениЯ формы
	V_x = 6;     // (м/с)
	V_y = 10;    // (м/с)
	m = 1;       // (кг)
	g = 9.81;    // (м/с^2)

	k = C*S*P/2; // безразмерная величина

	QScatterSeries* series1 = create_ser_without_resistance( V_x, V_y, g);
	QLineSeries* series2 = create_arbitrary_angle( m, k, V_x, V_y, m, g, 0, 0, 0, 0, 0);
	series2->setName("с учетом сопротивления(точное решение):");


	const int N = 20000;
	double _h = 0.001;

	QLineSeries* series3 = create_numer_approx ( N, _h, k, m, g, V_x, V_y);

	double V_x1 = 6,
		   V_y1 = 10;

	QLineSeries* series4 = create_relativism ( _h, V_y1, V_x1, m, k, g, N);

	chrt->addSeries(series1);
	chrt->addSeries(series2);
	chrt->addSeries(series3);
	chrt->addSeries(series4);

	QChartView *chartViw = new QChartView(this);
	chartViw->setMinimumSize(t.width() / 2, t.height());
	chartViw->setGeometry(t.width() / 2, 0, t.width() / 2, t.height());

	chart = new MyChart(0,25);

	chartViw->setChart(chart);
	chart->setTitle("Графики движения тела c непостоянной массой:");

	double m_p = 1,
		   m_p2 = 0.6,
		   V_xp = 6,  // начальная скорость
		   V_yp = 10,
		   V_xt = 0,
		   V_yt = 5,
		   dm = 0.000001,
		   t_detach = 0.6;

	QLineSeries* series5 = new QLineSeries();

	for (double t = 0.001; t < 20; t += 0.001) { //за промежуток времени dt откалывается кусочек массой dm по закону (в нуле откалывается 0)

		std::pair<double, double> coord = position_in_space( m_p, k, g, t, V_xp, V_yp);
		series5 = create_arbitrary_angle( m_p, k, V_xp, V_yp, dm, g, V_xt, V_yt, t, coord.first, coord.second);
		m_p = m_p - dm;

	}
	series5->setName("постепенное изменение массы:");

	QLineSeries* series6 = create_arbitrary_angle( m_p, k, V_xp, V_yp, m_p2, g, V_xt, V_yt, t_detach, 0, 0);
	series6->setName("случай резкого изменения массы под произвольным углом:");

	chart->addSeries(series5);
	chart->addSeries(series6);
}


MainWindow::~MainWindow()
{
	delete ui;
}


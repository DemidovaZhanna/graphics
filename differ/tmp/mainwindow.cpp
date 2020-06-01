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


void MainWindow::connect_axes( int a, int b, QChart *chrt, QXYSeries *series)
{
	axisX->setRange(a,b);
	axisX->setLabelFormat("%.2f");

	axisY->setRange(a,b);
	axisY->setLabelFormat("%g");

	chrt->addSeries(series);
	chrt->setAxisX(axisX, series);
	chrt->setAxisY(axisY, series);
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


QLineSeries* MainWindow::create_exact_solution ( double m, double V_x, double V_y, double k, double g)
{
	QLineSeries* series = new QLineSeries();
	series->setName("с учетом сопротивления(точное решение)");

	for (double t = 0; t < 20; t += 0.01)
		series->append( m*V_x/k*(1-exp(-k*t/m)), m/k*((V_y+m*g/k)*(1-exp(-k*t/m)) - g*t));

	return series;
}


QLineSeries* MainWindow::create_numer_approx ( int i, const int N, double _h, double k, double m, double g, double k1, double k2,
						  double k3, double k4, double l1, double l2, double l3, double l4, double pV_x[],
						  double pV_y[], double pY[], double pV[], double pX[])
{
	QLineSeries* series = new QLineSeries();
	series->setName("с учетом сопротивления(численное приближение)");

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

	return series;
}


QLineSeries* MainWindow::create_relativism( double k1, double l1, double k2, double l2, double k3, double l3, double k4,
					   double l4, double _h, double pV_y[], double pV_x[], double pX[], double pY[],
					   double pV[], double m, double k, double g, int j, const int N)
{
	QLineSeries* series = new QLineSeries();
	series->setName("на релятивистских скоростях(численное приближение)");

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

	return series;
}


QLineSeries* MainWindow::create_department_mass ( double x0, double y0, double m_p, double k, double dt, double V_xp,
									  double V_yp, double m_p2, double g, double V_t)
{
	QLineSeries* series = new QLineSeries();
	series->setName("случай резкого изменения массы");
	bool detached = false;

	for (double t = 0; t < 20; t += 0.01) {
		series->append(x0 + m_p*V_xp/k*(1-exp(-k*(t-dt)/m_p)), y0 + m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*(t-dt)/m_p)) - g*(t-dt)));

		if(t >= 0.6 && !detached){
			double V_x0 = V_xp*exp(-k/m_p * t);
			double V_y0 = V_yp*exp(-k/m_p * t) - g*m_p/k*(1 - exp(-k/m_p* t));
			x0 = m_p*V_xp/k*(1-exp(-k*t/m_p));
			y0 = m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*t/m_p)) - g*t);

			V_xp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_x0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
			V_yp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_y0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
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

	chrt = new QChart;
	axisX = new QValueAxis;
	axisY = new QValueAxis;

	chartView->setChart(chrt);
	chrt->setTitle("Графики движения тела под углом к горизонту:");

	connect_axes( 0, 13, chrt, series);

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

	QLineSeries* series2 = create_exact_solution ( m, V_x, V_y, k, g);


	const int N = 20000;
	double k1, k2, k3, k4, l1, l2, l3, l4;
	double pY[N]; pY[0] = 0;                                        // высота
	double pX[N]; pX[0] = 0;                                        //дальность
	double pV_x[N]; pV_x[0] = V_x;								    //массив - скорость: проекция на горизонталь
	double pV_y[N]; pV_y[0] = V_y;								    //массив - скорость: проекция на вертикаль
	double pV[N]; pV[0] = sqrt(pV_x[0]*pV_x[0] + pV_y[0]*pV_y[0]);  //полная скорость (величина вектора)
	int i = 1;
	double _h = 0.001;

	QLineSeries* series3 = create_numer_approx ( i, N, _h, k, m, g, k1, k2, k3, k4,
												 l1, l2, l3, l4, pV_x, pV_y, pY, pV, pX);

	double V_x1 = 6;
	double V_y1 = 10;																		
	double pV_x1[N]; pV_x1[0] = V_x1;										//массив - скорость: проекция на горизонталь
	double pV_y1[N]; pV_y1[0] = V_y1;										//массив - скорость: проекция на вертикаль
	double pV1[N]; pV1[0] = sqrt(pV_x1[0]*pV_x1[0] + pV_y1[0]*pV_y1[0]);    //полная скорость (величина вектора)

	QLineSeries* series4 = create_relativism( k1, l1, k2, l2, k3, l3, k4, l4, _h,
											  pV_y1, pV_x1, pX, pY, pV1, m, k, g, i, N);

	chrt->addSeries(series1);
	chrt->addSeries(series2);
	chrt->addSeries(series3);
	chrt->addSeries(series4);

	chrt->setAxisX(axisX, series1);
	chrt->setAxisY(axisY, series1);

	chrt->setAxisX(axisX, series2);
	chrt->setAxisY(axisY, series2);

	chrt->setAxisX(axisX, series3);
	chrt->setAxisY(axisY, series3);

	chrt->setAxisX(axisX, series4);
	chrt->setAxisY(axisY, series4);


	QChartView *chartViw = new QChartView(this);
	chartViw->setMinimumSize(t.width() / 2, t.height());
	chartViw->setGeometry(t.width() / 2, 0, t.width() / 2, t.height());
	chart = new QChart;

	chartViw->setChart(chart);
	chart->setTitle("Графики движения тела c непостоянной массой:");

	QValueAxis *axisX1 = new QValueAxis;
	axisX1->setRange(0,25 );
	axisX1->setLabelFormat("%.2f");

	QValueAxis *axisY1 = new QValueAxis;
	axisY1->setRange(0,13);
	axisY1->setLabelFormat("%g");

	double m_p = 1;
	double m_p2 = 0.6;
	double m_t = 0.1;
	double V_xp = 6; // начальная скорость
	double V_yp = 10;
	double V_t = 0.5; //скорость топлива
	QLineSeries* series5 = new QLineSeries();
	series5->setName("постепенное изменение массы");

	for (double t = 0; t < 20; t += 0.01) {
//		double m1 = m_p - m_t/(1 + t*t);
//		series5->append( m1*V_xp/k*(1-exp(-k*t/m1)), m1/k*((V_yp+m*g/k)*(1-exp(-k*t/m1)) - g*t));
	}

	QLineSeries* series6 = new QLineSeries();
	series6->setName("случай резкого изменения массы");
	bool detached = false;
	double dt = 0;
	double x0 = 0, y0 = 0;

	for (double t = 0; t < 20; t += 0.01) {
		series6->append(x0 + m_p*V_xp/k*(1-exp(-k*(t-dt)/m_p)), y0 + m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*(t-dt)/m_p)) - g*(t-dt)));

		if(t >= 0.6 && !detached){
			double V_x0 = V_xp*exp(-k/m_p * t);
			double V_y0 = V_yp*exp(-k/m_p * t) - g*m_p/k*(1 - exp(-k/m_p* t));
			x0 = m_p*V_xp/k*(1-exp(-k*t/m_p));
			y0 = m_p/k*((V_yp+m_p*g/k)*(1-exp(-k*t/m_p)) - g*t);

			V_xp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_x0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
			V_yp = ((m_p*sqrt(V_x0*V_x0 + V_y0*V_y0) + (m_p - m_p2)*V_t)/m_p2)*V_y0*sqrt(1/(V_y0*V_y0 + V_x0*V_x0));
			m_p = m_p2;
			detached = true;
			dt = t;
		}
	}

	chart->addSeries(series5);
	chart->addSeries(series6);

	chart->setAxisX(axisX1, series5);
	chart->setAxisY(axisY1, series5);

	chart->setAxisX(axisX1, series6);
	chart->setAxisY(axisY1, series6);


}


MainWindow::~MainWindow()
{
	delete ui;
}


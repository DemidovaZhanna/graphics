#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QDesktopWidget>
#include <math.h>

double f_V_x(double V_x, double k, double m)
{
	return -k/m*V_x;
}

double f_V_y(double V_y, double k, double m, double g)
{
	return -k/m*V_y- g;
}

double f_Vx1(double V_x1, double V_y1, double k, double m)
{
	return -k*V_x1*sqrt(V_x1*V_x1 + V_y1*V_y1)/m;
}

double f_Vy1(double V_x1, double V_y1, double k, double m, double g)
{
	return -k*V_y1*sqrt(V_x1*V_x1 + V_y1*V_y1)/m - g;
}

void connect_axes(int a, int b, QChart *chrt, QXYSeries *series)
{
	chrt->addSeries(series);
	QValueAxis *axisX = new QValueAxis;
	axisX->setRange(a,b);
	axisX->setLabelFormat("%.2f");

	QValueAxis *axisY = new QValueAxis;
	axisY->setRange(a,b);
	axisY->setLabelFormat("%g");

	chrt->setAxisX(axisX, series);
	chrt->setAxisY(axisY, series);
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
	QValueAxis *axisX = new QValueAxis;
	axisX->setRange(0,13);
	axisX->setLabelFormat("%.2f");

	QValueAxis *axisY = new QValueAxis;
	axisY->setRange(0,13);
	axisY->setLabelFormat("%g");

	chartView->setChart(chrt);
	chrt->setTitle("Графики движения тела под углом к горизонту:");

	double m, V_x, V_y, C, S, P, k, g;

	P = 1.29;    // (кг/м^3)
	S = 0.1;     // S = M_PI*r^2
	C = 0.47;    // безразмерный коэффициент сопротивлениЯ формы
	V_x = 6;     // (м/с)
	V_y = 10;    // (м/с)
	m = 1;       // (кг)
	g = 9.81;    // (м/с^2)

	k = C*S*P/2; // безразмерная величина

	QScatterSeries* series1 = new QScatterSeries();
	series1->setMarkerSize(10);
	series1->setName("без учета сопротивления воздуха");

	QLineSeries* series2 = new QLineSeries();
	series2->setName("с учетом сопротивления(точное решение)");

	QLineSeries* series3 = new QLineSeries();
	series3->setName("с учетом сопротивления(численное приближение)");

	QLineSeries* series4 = new QLineSeries();
	series4->setName("на релятивистских скоростях(численное приближение)");

	for (double t = 0; t < 20; t += 0.01) {
		series1->append( V_x*t,V_y*t-g*t*t/2);
		series2->append( m*V_x/k*(1-exp(-k*t/m)), m/k*((V_y+m*g/k)*(1-exp(-k*t/m)) - g*t));
	}

	const int N = 20000;
	double k1, k2, k3, k4, l1, l2, l3, l4;
	double pY[N]; pY[0]=0; // высота
	double pX[N]; pX[0]=0; //дальность
	//double pt[N]; pt[0]=0; //время
	double pV_x[N]; pV_x[0]=V_x; //массив - скорость: проекция на горизонталь
	double pV_y[N]; pV_y[0]=V_y; //массив - скорость: проекция на вертикаль
	double pV[N]; pV[0]= sqrt(pV_x[0]*pV_x[0] + pV_y[0]*pV_y[0]); //полная скорость (величина вектора)
	int i = 1;
	double _h = 0.001;

	while(pY[i-1]>-0.00001 && i < N)
		{
			//pt[i] =pt[i-1]+_h;
			k1 = _h * f_V_x( pV_x[i-1],k, m);
			l1 = _h * f_V_y( pV_y[i-1], k, m, g);
			k2 = _h * f_V_x( pV_x[i-1]+k1/2, k, m);
			l2 = _h * f_V_y( pV_y[i-1]+l1/2, k, m, g);
			k3 = _h * f_V_x( pV_x[i-1]+k2/2, k, m);
			l3 = _h * f_V_y( pV_y[i-1]+l2/2, k, m, g);
			k4 = _h * f_V_x( pV_x[i-1]+k3, k, m);
			l4 = _h * f_V_y( pV_y[i-1]+l3, k, m, g);


			// получаем проекции скорости из предварительно полученных приближенных значений
			pV_x[i]=pV_x[i-1] + ( k1+2*k2+2*k3+k4)/6;
			pV_y[i]=pV_y[i-1] + ( l1+2*l2+2*l3+l4)/6;

			pV[i] = sqrt( pV_x[i]*pV_x[i] + pV_y[i]*pV_y[i] ); //величина вектора скорости из проекций

			pY[i]=pY[i-1]+pV_y[i]*_h;
			pX[i]=pX[i-1]+pV_x[i]*_h; //пересчет координат

			series3->append(pX[i], pY[i]);

			i++;
		}

	double V_x1 = 6;
	double V_y1 = 10;
	double m1, m2, m3, m4, f1, f2, f3, f4;
	int j =1;

	double pY1[N]; pY1[0]=0; // высота
	double pX1[N]; pX1[0]=0; //дальность
	//double pt1[N]; pt1[0]=0; //время
	double pV_x1[N]; pV_x1[0]=V_x1; //массив - скорость: проекция на горизонталь
	double pV_y1[N]; pV_y1[0]=V_y1; //массив - скорость: проекция на вертикаль
	double pV1[N]; pV1[0]= sqrt(pV_x1[0]*pV_x1[0] + pV_y1[0]*pV_y1[0]); //полная скорость (величина вектора)

	while (pY1[j-1]>-0.00001 && j < N )
		{
			//pt1[i] = pt1[i-1]+_h;
			m1 = _h * f_Vx1( pV_y1[j-1], pV_x1[j-1], k, m);
			f1 = _h * f_Vy1( pV_y1[j-1], pV_x1[j-1], k, m, g);
			m2 = _h * f_Vx1( pV_y1[j-1]+f1/2, pV_x1[j-1]+m1/2, k, m);
			f2 = _h * f_Vy1( pV_y1[j-1]+f1/2, pV_x1[j-1]+m1/2, k, m, g);
			m3 = _h * f_Vx1( pV_y1[j-1]+f2/2, pV_x1[j-1]+m2/2, k, m);
			f3 = _h * f_Vy1( pV_y1[j-1]+f2/2, pV_x1[j-1]+m2/2, k, m, g);
			m4 = _h * f_Vx1( pV_y1[j-1]+f3, pV_x1[j-1]+m3, k, m);
			f4 = _h * f_Vy1( pV_y1[j-1]+f3, pV_x1[j-1]+m3, k, m, g);

			// получаем проекции скорости из предварительно полученных приближенных значений
			pV_x1[j]=pV_x1[j-1] + (m1+2*m2+2*m3+m4)/6;
			pV_y1[j]=pV_y1[j-1] + (f1+2*f2+2*f3+f4)/6;

			pV1[j] = sqrt( pV_x1[j]*pV_x1[j] + pV_y1[j]*pV_y1[j]); //величина вектора скорости из проекций

			pY1[j]=pY1[j-1]+pV_y1[j]*_h;
			pX1[j]=pX1[j-1]+pV_x1[j]*_h; //пересчет координат

			series4->append(pX1[j], pY1[j]);

			j++;
		}

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
	axisX->setLabelFormat("%.2f");

	QValueAxis *axisY1 = new QValueAxis;
	axisY1->setRange(0,13);
	axisY->setLabelFormat("%g");

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


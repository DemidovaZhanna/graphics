#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include "my_chart.h"

QT_BEGIN_NAMESPACE
using namespace QtCharts;

namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();
	QLineSeries* create_exact_solution ( double m, double V_x, double V_y, double k, double g, double x0, double y0);
	QScatterSeries* create_ser_without_resistance(double V_x, double V_y, double g, double x0, double y0);
	QLineSeries* create_numer_approx (const int N, double t, double k, double m, double g, double V_x,
									  double V_y, double x0, double y0);
	QLineSeries* create_relativism(double _h, double V_y, double V_x, double m, double k, double g, int N, double x0, double y0);
	QLineSeries* create_slow_change (double m_p, double k, double V_xp,
												 double V_yp, double m_t, double g, double V_xt, double V_yt);

//	QLineSeries* create_detaching_mass (double m_p, double k, double V_xp,
//										 double V_yp, double m_p2, double g, double V_t);
	QLineSeries* create_arbitrary_angle (double m_p, double k, double V_xp,
										 double V_yp, double m_p2, double g, double V_xt, double V_yt,
										 double t, double x_det, double y_det, double x0, double y0);


private:
	Ui::MainWindow *ui;

	MyChart *chrt;
	MyChart *chart;

	QValueAxis *axisX, *axisY;
};

#endif // MAINWINDOW_H

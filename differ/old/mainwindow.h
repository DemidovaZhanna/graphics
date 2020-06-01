#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>

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
	void connect_axes(int a, int b);

private:
	Ui::MainWindow *ui;

	QChart *chrt;
	QChart *chart;
};
#endif // MAINWINDOW_H

#include "Database.h"
#include "Helper_Functions.h"
#include <gtest/gtest.h>

namespace {

using ::testing::TestWithParam;
using ::testing::Values;

TEST(DatabaseTest, Vectors) {
    string config_file = "./conf/cluster.conf";
    string input_file =  "./datasets/test_database_vectors.csv";
    /* read test data */
    int* cluster_config = new int[4];
    vector<vector<double>> cluster_data;
    vector<string> item_ids;
    ASSERT_EQ(1, Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids));
    /* build database */
    DistanceDatabase<double>* db = new DistanceDatabase<double>();
    db->calculate_distances(&cluster_data);
    /* run queries */
    EXPECT_DOUBLE_EQ(301.77927936333538, db->get_distance(0, 1));
    EXPECT_DOUBLE_EQ(301.77927936333538, db->get_distance(1, 0));
    /* free memory */
    delete (db);
}

TEST(DatabaseTest, Curves) {
    string config_file = "./conf/cluster.conf";
    string input_file =  "./datasets/test_database_curves.csv";
    /* read test data */
    int* cluster_config = new int[4];
    vector<vector<double*>> cluster_data;
    vector<string> item_ids;
    ASSERT_EQ(1, Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids));
    /* build database */
    DistanceDatabase<double*>* db = new DistanceDatabase<double*>();
    db->calculate_distances(&cluster_data);
    /* run queries */
    EXPECT_DOUBLE_EQ(0.57986771559476025, db->get_distance(0, 1));
    EXPECT_DOUBLE_EQ(0.57986771559476025, db->get_distance(1, 0));
    /* free memory */
    delete (db);
}

TEST(DatabasTest, OutOfBounds) {
    string config_file = "./conf/cluster.conf";
    string input_file =  "./datasets/test_database_curves.csv";
    /* read test data */
    int* cluster_config = new int[4];
    vector<vector<double*>> cluster_data;
    vector<string> item_ids;
    ASSERT_EQ(1, Read_files(&cluster_data, cluster_config, input_file, config_file, &item_ids));
    /* build database */
    DistanceDatabase<double*>* db = new DistanceDatabase<double*>();
    db->calculate_distances(&cluster_data);
    /* run queries */
    EXPECT_DOUBLE_EQ(-1.0, db->get_distance(-1, 1));
    EXPECT_DOUBLE_EQ(-1.0, db->get_distance(-100, -4));
    EXPECT_DOUBLE_EQ(-1.0, db->get_distance(5, 0));
    EXPECT_DOUBLE_EQ(-1.0, db->get_distance(721, 256));
    EXPECT_DOUBLE_EQ(-1.0, db->get_distance(56, 56));
    /* free memory */
    delete (db);
}

}


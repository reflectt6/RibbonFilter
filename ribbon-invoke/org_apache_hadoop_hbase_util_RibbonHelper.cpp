//
// Created by Rain Night on 2024/1/11.
//
#include "org_apache_hadoop_hbase_util_RibbonHelper.h"
#include <chrono>

using TypeParam = tmp_Settings;
IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
IMPORT_RIBBON_IMPL_TYPES(TypeParam);

// 全局变量声明
InterleavedSoln *isoln = nullptr;
Banding *banding = nullptr;
Hasher *hasher = nullptr;
double addDurationMs;
double queryDurationMs;
double backSubstMs;
double stringToCharsMs;
double initMs;

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    initRibbonFilter
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_initRibbonFilter
        (JNIEnv *, jobject, jint ji) {
    auto start_time = std::chrono::high_resolution_clock::now();

//    Index num_slots = 131072;
    Index num_slots = static_cast<uint32_t>(ji);
//    uint32_t num_to_add = 109306;
    uint32_t max_ibytes = static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
    std::unique_ptr<char[]> idata(new char[max_ibytes]);
    isoln = new InterleavedSoln(idata.get(), max_ibytes);
    banding = new Banding();
    banding->SetOrdinalSeed(0);
    banding->Reset(num_slots);
    hasher = new Hasher();
    addDurationMs = 0.0;
    queryDurationMs = 0.0;
    backSubstMs = 0.0;
    stringToCharsMs = 0.0;
    initMs = 0.0;


    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double duration_ms = static_cast<double>(duration.count()) / 1000.0;
    initMs += duration_ms;
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    addKey
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_addKey
        (JNIEnv * env, jobject jo, jstring js) {
    auto start_time = std::chrono::high_resolution_clock::now();
    const char *cString = env->GetStringUTFChars(js , nullptr);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double duration_ms = static_cast<double>(duration.count()) / 1000.0;
    stringToCharsMs += duration_ms;

    start_time = std::chrono::high_resolution_clock::now();
    bool b = banding->Add(Slice(cString));
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    duration_ms = static_cast<double>(duration.count()) / 1000.0;
    addDurationMs += duration_ms;
    return b;
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    backSubst
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_backSubst
        (JNIEnv *, jobject) {
    auto start_time = std::chrono::high_resolution_clock::now();
    isoln->BackSubstFrom(*banding);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double duration_ms = static_cast<double>(duration.count()) / 1000.0;
    backSubstMs += duration_ms;
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    filterQuery
 * Signature: (Ljava/lang/String;)Z
 */
JNIEXPORT jboolean JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_filterQuery
        (JNIEnv * env, jobject jo, jstring js) {
    auto start_time = std::chrono::high_resolution_clock::now();
    const char *cString = env->GetStringUTFChars(js , nullptr);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double duration_ms = static_cast<double>(duration.count()) / 1000.0;
    stringToCharsMs += duration_ms;


    start_time = std::chrono::high_resolution_clock::now();
    bool b = isoln->FilterQuery(Slice(cString), *hasher);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    duration_ms = static_cast<double>(duration.count()) / 1000.0;
    queryDurationMs += duration_ms;
    return b;
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    printDuration
 * Signature: ()Z
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_printDuration
        (JNIEnv *, jobject) {
    printf("ribbon native add Key time = %.2f ms\n",addDurationMs);
    printf("ribbon native backSubst time = %.2f ms\n",backSubstMs);
    printf("ribbon native query time = %.2f ms\n",queryDurationMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    getAddDuration
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_getAddDuration
        (JNIEnv *, jobject) {
    return static_cast<jdouble>(addDurationMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    getBackSubstDuration
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_getBackSubstDuration
        (JNIEnv *, jobject) {
    return static_cast<jdouble>(backSubstMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    getQueryDuration
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_getQueryDuration
        (JNIEnv *, jobject) {
    return static_cast<jdouble>(queryDurationMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    getStringToCharsDuration
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_getStringToCharsDuration
        (JNIEnv *, jobject) {
    return static_cast<jdouble>(stringToCharsMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    getInitDuration
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_getInitDuration
        (JNIEnv *, jobject) {
    return static_cast<jdouble>(initMs);
}

/*
 * Class:     org_apache_hadoop_hbase_util_RibbonHelper
 * Method:    close
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_apache_hadoop_hbase_util_RibbonHelper_close
        (JNIEnv *, jobject) {
    delete banding;
    delete isoln;
    delete hasher;
}